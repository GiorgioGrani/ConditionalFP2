import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.util.List;
import java.util.TreeMap;
import java.util.Map;

import java.util.ArrayList;

public class RFP  {
    private double[] c;
    private double[][] A;
    private double [] b;
    private double [] lower;
    private double [] upper;
    private int numberOfIntegerVariables;
    private int n ;
    private static final double eps = 1e-6;

    private IloCplex cplex;
    private IloModel model;

    private IloAddable objective;
    private IloAddable reductionobj;
    private ArrayList< IloNumVar> x;
    private ArrayList<IloNumVar> r;
    private ArrayList<IloNumVar> e;
    private IloLinearNumExpr halfobj;
    private IloLinearNumExpr [] leftcons;
    private IloAddable[] constraints;
    private double constant;
    private IloAddable distconstraint;
    private int count;
    private double [] v;

    public void set(double[] c, double[][] A , double [] b, double [] lower , double [] upper, int numberOfIntegerVariables) throws IloException{
        if( c != null) this.c = c;
        else return;
        this.n = this.c.length;
        if( A != null){
            this.A = A;
            if( c.length != A[0].length) return;
        }

        if( b != null) this.b = b;

        if( lower != null)
            this.lower = lower;
        else
            this.lower = this.fill(false);

        if( upper != null)
            this.upper = upper;
        else
            this.upper = this.fill(true);

        this.numberOfIntegerVariables = numberOfIntegerVariables;

        this.cplex = new IloCplex();
        this.model = this.cplex.getModel();

        this.createVariables();
        this.setObjective();
        this.setConstraints();
    }

    private void createVariables() throws IloException {
        this.x = new ArrayList<>();
        for ( int i = 0 ; i< this.n; i++){
            String varname = "x" + i;
            IloNumVar xi = this.cplex.numVar(this.lower[i], this.upper[i], varname);
            this.x.add( xi);
        }
    }


    //todo rimuovere objective dopo averlo utilizzato



    private void setObjective() throws IloException{

        IloLinearNumExpr expr = this.cplex.linearNumExpr();
        int i = 0;
        for(IloNumVar xi : x){
            expr.addTerm(xi, this.c[i]);
            i++;
        }
        this.objective = this.cplex.addMinimize(expr);
    }

    private void setConstraints() throws IloException{
        if(this.A == null || this.b == null) return;

        IloAddable[] constraints = new IloAddable[this.A.length];

        this.leftcons = new IloLinearNumExpr[this.A.length];
        for( int i = 0 ; i < this.A.length ; i++){
            int j = 0;
            IloLinearNumExpr expr = this.cplex.linearNumExpr();
            for( IloNumVar xi : x){
                expr.addTerm( A[i][j], xi);
                j++;
            }
            this.leftcons[i] = expr;
            constraints[i] = (this.cplex.addGe(expr, b[i]));
        }
        this.constraints = constraints;
    }




    private double[] fill(boolean up){
        double [] d = new double [this.n];
        for(int i = 0 ; i < n; i++){
            d[i] = (up? Double.MAX_VALUE : -(Double.MAX_VALUE -1));
        }
        return d;
    }

    private ArrayList<Object> getXTilde() throws IloException{
        ArrayList<Double> ret = new ArrayList<>();
        ArrayList<Double> x = new ArrayList<>();
        boolean stop = true;
        for(int i = 0; i < this.numberOfIntegerVariables; i++){
            IloNumVar xi = this.x.get(i);
            Double d = this.cplex.getValue(xi);
            Double rd = FButils.round(d);
            ret.add(rd);
            x.add(d);
            if(Math.abs(rd - d) >= this.eps){
                stop = false;
            }
        }
        ArrayList<Object> res = new ArrayList<>();
        res.add(ret);
        res.add(stop);
        res.add(x);
        return res;
    }
    private ArrayList<Double> getX() throws IloException{
        ArrayList<Double> x = new ArrayList<>();
        for(int i = 0; i < this.numberOfIntegerVariables; i++){
            Double d = this.cplex.getValue(this.x.get(i));
            x.add(d);
        }
        return x;
    }


    private boolean checkFeasibility(ArrayList<Double> x){
        if(this.A == null || this.b == null) return true;
        for(int i = 0 ; i < this.A.length; i ++){
            double val = 0;
            for( int j = 0 ; j < this.A[0].length; j++){
                double xi = x.get(j);
                val += xi*A[i][j];
            }
            if (val < b[i] - eps) return false;
        }
        return true;
    }

    private double checkGap(ArrayList<Double> x, double lpobj){
        double val = 0;
        for ( int i = 0 ; i < this.c.length; i++){
            val += c[i]*x.get(i);
        }
        return - lpobj + val;
    }

    private boolean stopCondition(ArrayList<Double> x) {
        return FButils.integralityGap(IntegralityGapTypes.L2Norm,this.numberOfIntegerVariables, x) <= FButils.eps;
    }

    private void createDistanceVariables() throws IloException{
        ArrayList<IloNumVar> e = new ArrayList<>();
        ArrayList<IloNumVar> r = new ArrayList<>();
        for(int i = 0; i < this.A.length; i++){
            e.add(this.cplex.numVar(0, Double.MAX_VALUE, "e"+i));
            r.add(this.cplex.numVar(0, Double.MAX_VALUE, "r"+i));
        }
        this.e = e;
        this.r = r;
    }

    private void buildHalfRedObj() throws IloException{
        IloLinearNumExpr expr = this.cplex.linearNumExpr();
        for (int i = 0; i< this.A.length; i++){
            expr.addTerm(-1d, this.e.get(i));
            expr.addTerm(-1d, this.r.get(i));
        }
        this.halfobj = expr;
    }

    private void setReductionConstraints() throws IloException {
        this.model.remove(this.constraints);
        for (int i = 0; i < this.A.length; i++) {
            IloLinearNumExpr expr = this.cplex.linearNumExpr();
            expr.addTerm(1d, this.e.get(i));
            IloNumExpr sum = this.cplex.sum(expr, b[i]);
            this.cplex.addGe(this.leftcons[i], sum);
        }

        this.v = this.setV();

        for (int i = 0; i < this.A.length; i++) {
            IloLinearNumExpr expr = this.cplex.linearNumExpr();
            expr.addTerm(-1d, this.r.get(i));
            IloNumExpr sum = this.cplex.sum(expr, v[i]);
            this.cplex.addLe(this.leftcons[i], sum);
        }

    }


    private double[] setV() throws IloException{
        IloCplex subcplex = new IloCplex();
        IloModel submodel = subcplex.getModel();
        subcplex.setOut(null);

        ///// variables
        ArrayList<IloNumVar> x = new ArrayList<>();
        for(int i = 0; i < this.n; i++){
            x.add(subcplex.numVar(this.lower[i], this.upper[i]));
        }

        ///// constraints
        ArrayList<IloLinearNumExpr> objs = new ArrayList<>();
        for(int i = 0; i < this.A.length; i++){
            IloLinearNumExpr expr = subcplex.linearNumExpr();
            for(int j = 0; j < this.A[0].length; j++){
                expr.addTerm(A[i][j], x.get(i));
            }
            objs.add(expr);
            subcplex.addGe(expr,this.b[i]);
        }

        ////// finding v_i
        double [] v = new double [this.A.length];
        double constant = 0;
        int i = 0;
        for( IloLinearNumExpr obj : objs){
            IloAddable currentobj = subcplex.addMaximize(obj);
            subcplex.solve();
            v[i] = subcplex.getObjValue();
            constant += v[i] - b[i];
            i++;
            submodel.remove(currentobj);
        }
        this.constant = constant;

        return v;
    }

    private void setReductionObjective(ArrayList<Double> x) throws IloException{
        IloLinearNumExpr expr = this.cplex.linearNumExpr();
        int i = 0;
        int count = 0;
        for( Double xi : x){
            if(xi <= 0.5) expr.addTerm(this.constant, this.x.get(i));
            else{
                count ++;
                expr.addTerm(-this.constant, this.x.get(i));
            }
            i++;
        }
        this.count = count;

        this.model.remove(this.distconstraint);
        this.distconstraint = this.cplex.addGe(expr, this.constant*(1d - count));
        this.model.remove(this.reductionobj);
        IloNumExpr obj = this.cplex.sum( this.halfobj, expr);
        this.reductionobj = this.cplex.addMinimize(obj);
    }
    public ArrayList<Object> solve() throws IloException{
        int maxiter = 0;
        this.cplex.setOut(null);

        long start = System.currentTimeMillis();
        this.cplex.solve();
        ArrayList<Object> ret = new ArrayList<>();
        ArrayList<Double> x = this.getX();
        boolean stop = this.stopCondition(x);
        double lpobj = this.cplex.getObjValue();
        if(stop){
            long end = System.currentTimeMillis();
            ret.add(x);
            ret.add(stop);
            ret.add(0);
            ret.add(-start + end);
            ret.add(this.checkFeasibility(x));
            ret.add( this.checkGap(x,lpobj));
            return ret;
        }

        this.model.remove(this.objective);
        this.createDistanceVariables();
        this.setReductionConstraints();
        this.buildHalfRedObj();
        while(!stop && maxiter <= 100) {
            maxiter++;
            this.setReductionObjective(x);
            this.cplex.solve();
            x = this.getX();
            stop = this.stopCondition(x);
            double oval = this.cplex.getObjValue();
            double hval = this.cplex.getValue(this.halfobj);
            double dist = oval - hval;
            //System.out.println(this.constant+"  "+maxiter+") "+FButils.integralityGap(IntegralityGapTypes.L2Norm,this.numberOfIntegerVariables,x)+"  dist:"+dist+
              //      "  normdist:"+(count+(dist/this.constant))+"  hval: "+ hval+"   tot: "+ oval);
//            for(int i = 0; i< this.A.length; i++){
//                System.out.println("              "+(this.v[i] - this.cplex.getValue(this.r.get(i)))+"  "+ (this.b[i] + this.cplex.getValue(this.e.get(i))));//this.model.remove(this.distanceobj);
//            }
        }
        long end = System.currentTimeMillis();
        ret.add(x);
        ret.add(stop);
        ret.add(maxiter-1);
        ret.add(-(start - end));
        ret.add(this.checkFeasibility(x));
        ret.add( this.checkGap(x,lpobj));
        return ret;
    }

}
