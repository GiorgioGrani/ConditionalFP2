import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.util.List;
import java.util.TreeMap;
import java.util.Map;

import java.util.ArrayList;

public class Simple  {
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
            if(xi <= 0.5) expr.addTerm(1d, this.x.get(i));
            else expr.addTerm(-1d, this.x.get(i));
            i++;
            if( i >= this.numberOfIntegerVariables) break;
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
    private boolean stopCondition(ArrayList<Double> x, ArrayList<Double> xold){
        double val = 0;
        boolean equals = false;
        int i = 0;
        for(Double xi : x){
            if(Math.abs(xi - xold.get(i)) >= FButils.eps){
                equals = false;
                break;
            }else {
                equals = true;
            }
            i++;
        }
        return equals || FButils.integralityGap(IntegralityGapTypes.L2Norm,this.numberOfIntegerVariables, x) <= FButils.eps;
    }

    public ArrayList<Object> solve() throws IloException{
        int maxiter = 0;
        this.cplex.setOut(null);

        long start = System.currentTimeMillis();
        this.cplex.solve();
        ArrayList<Object> ret = new ArrayList<>();
        ArrayList<Double> x = this.getX();
        ArrayList<Double> xold = x;
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
        while(!stop && maxiter <= 100) {
            maxiter++;
            this.setDistanceObjective(x);
            this.cplex.solve();
            x = this.getX();
            stop = this.stopCondition(x, xold);
            //System.out.println(maxiter+") "+FButils.integralityGap(IntegralityGapTypes.L2Norm,this.numberOfIntegerVariables,x));
            this.model.remove(this.distanceobj);
            xold = x;
        }
        long end = System.currentTimeMillis();
        ret.add(x);
        ret.add(stop);
        ret.add(maxiter-1);
        ret.add(-(start - end));
        ret.add(this.checkFeasibility(x));
        ret.add( FButils.objVal(x, this.c));
        return ret;
    }

}
