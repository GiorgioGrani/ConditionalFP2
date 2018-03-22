import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.util.List;
import java.util.TreeMap;
import java.util.Map;

import java.util.ArrayList;

public class LJ4S  {
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
    private IloAddable distanceobj;
    private ArrayList< IloNumVar> x;

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
            IloNumVar xi = this.cplex.numVar(this.lower[i], this.upper[i]);
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

        for( int i = 0 ; i < this.A.length ; i++){
            int j = 0;
            IloLinearNumExpr expr = this.cplex.linearNumExpr();
            for( IloNumVar xi : x){
                expr.addTerm( A[i][j], xi);
                j++;
            }
            this.cplex.addGe(expr, b[i]);
        }
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

    private ArrayList<Double> getXTilde(ArrayList<Double> x) throws IloException{
        ArrayList<Double> ret = new ArrayList<>();
        boolean stop = true;
        for(int i = 0; i < this.numberOfIntegerVariables; i++){
            Double xi = x.get(i);
            Double rd = FButils.round(xi);
            ret.add(rd);
            //System.out.println(rd+"   "+xi);
        }
        return ret;
    }
    private ArrayList<Double> getX() throws IloException{
        ArrayList<Double> x = new ArrayList<>();
        for(int i = 0; i < this.numberOfIntegerVariables; i++){
            Double d = this.cplex.getValue(this.x.get(i));
            x.add(d);
        }
        return x;
    }

    private void setDistanceObjective(ArrayList<Double> x) throws IloException{
        IloLinearNumExpr expr = this.cplex.linearNumExpr();
        int i = 0;
        for( Double xi : x){
            double d = -1.5*Math.PI*Math.cos(1.5*Math.PI*(xi + 1d));
            expr.addTerm(d, this.x.get(i));
            i++;
            if(i >= this.numberOfIntegerVariables) break;
        }
        this.distanceobj = this.cplex.addMinimize(expr);
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

    private boolean stopCondition(ArrayList<Double> x){
        return FButils.integralityGap(IntegralityGapTypes.L2Norm,this.numberOfIntegerVariables, x) <= FButils.eps;
    }

    private ArrayList<Double> getDirection(ArrayList<Double> x, ArrayList<Double> xk){
        ArrayList<Double> dk = new ArrayList<>();
        int i = 0;
        for(Double xi : x){
            if(i >= this.numberOfIntegerVariables) break;
            double val = xi- xk.get(i);
            dk.add(val);
            i++;
        }
        return dk;
    }

    private ArrayList<Double> nextStep(ArrayList<Double> xk, double ak, ArrayList<Double> dk){
        ArrayList<Double> step = new ArrayList<>();
        int i = 0;
        for(Double xi : xk){
            if(i >= this.numberOfIntegerVariables) break;
            double val = xi + ak*dk.get(i);
            step.add(val);
            i++;
        }
        return step;
    }

    private double Armijo(ArrayList<Double> xk, ArrayList<Double> dk){
        double ak = 1;
        double g = 1e-4;
        double d = 0.9;

        ArrayList<Double> grad = grad(xk);
        double value = scalarProd(dk, grad);
        while( f(nextStep(xk,ak,dk)) > f(xk) + g*ak*value){
            ak = ak*d;
        }
        return ak;
    }
    private double Armijo(ArrayList<Double> xk, ArrayList<Double> dk, double delta){
        double ak = delta;
        double g = 1e-4;
        double d = 0.9;

        ArrayList<Double> grad = grad(xk);
        double value = scalarProd(dk, grad);
        while( f(nextStep(xk,ak,dk)) > f(xk) + g*ak*value){
            ak = ak*d;
        }
        return ak;
    }

    private double scalarProd(ArrayList<Double> a, ArrayList<Double> b){
        double ret = 0;
        int i = 0;
        for(Double ai : a){
            ret = ret + ai*b.get(i);
            i++;
        }
        return ret;
    }

    private double f(ArrayList<Double> x){
        double ret = 0;
        for(Double xi : x){
            ret = ret + Math.sin(1.5*Math.PI*(xi + 1d))+1;
        }
        return ret;
    }

    private ArrayList<Double> grad(ArrayList<Double> x){
        ArrayList<Double> ret = new ArrayList<Double>();
        for(Double xi : x){
            double val = -1.5*Math.PI*Math.cos(1.5*Math.PI*(xi + 1d));;
            ret.add(val);
        }
        return ret;
    }

    public ArrayList<Object> solve() throws IloException{
        int maxiter = 0;
        this.cplex.setOut(null);

        long start = System.currentTimeMillis();
        this.cplex.solve();
        ArrayList<Object> ret = new ArrayList<>();
        ArrayList<Double> x = this.getX();
        ArrayList<Double> xk = x;
        boolean stop = this.stopCondition(xk);
        boolean round = false;
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
        while(!stop && maxiter <= 100) {
            maxiter++;
            this.setDistanceObjective(xk);
            this.cplex.solve();
            double val = this.cplex.getObjValue() - scalarProd(grad(xk), xk);
            if( -val <= 1e-3){ //todo: hai fissato questo stop molto lasco
                stop = true;
                round = true;
                break;
            }
            x = this.getX();
            ArrayList<Double> dk = this.getDirection(x, xk);
            double ak = 1d/(maxiter+0d);//
            //double ak =  this.Armijo(xk, dk);
            xk = this.nextStep(xk, ak , dk);
            stop = this.stopCondition(xk);
            //System.out.println(maxiter+") "+val +"    "+FButils.integralityGap(IntegralityGapTypes.L2Norm,this.numberOfIntegerVariables,xk));
            this.model.remove(this.distanceobj);
        }
        if(!round) {
            long end = System.currentTimeMillis();
            ret.add(xk);
            ret.add(stop);
            ret.add(maxiter - 1);
            ret.add(-(start - end));
            ret.add(this.checkFeasibility(xk));
            ret.add(FButils.objVal(xk, this.c));
        }else{
            xk = this.getXTilde(x);
            long end = System.currentTimeMillis();
            ret.add(xk);
            ret.add(stop);
            ret.add(maxiter - 1);
            ret.add(-(start - end));
            ret.add(this.checkFeasibility(xk));
            ret.add(FButils.objVal(xk, this.c));
        }
        return ret;
    }

}
