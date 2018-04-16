import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.util.List;
import java.util.TreeMap;
import java.util.Map;

import java.util.ArrayList;

public class BenchMark  {
    private double[] c;
    private double[][] A;
    private double [] b;
    private double [] lower;
    private double [] upper;
    private int [] directions;
    private int numberOfIntegerVariables;
    private int n ;
    private static final double eps = 1e-6;

    private IloCplex cplex;
    private IloModel model;

    private IloAddable objective;
    private IloAddable distanceobj;
    private ArrayList< IloNumVar> x;

    public void set(double[] c, double[][] A , double [] b, double [] lower , double [] upper, int numberOfIntegerVariables, int [] directions) throws IloException{
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

        this.directions = directions;
        this.cplex = new IloCplex();
        this.model = this.cplex.getModel();

        this.createVariables();
        this.setObjective();
        this.setConstraints();
    }

    private void createVariables() throws IloException {
        this.x = new ArrayList<>();
        for ( int i = 0 ; i< this.numberOfIntegerVariables; i++){
            String varname = "x" + i;
            int low = (int) Math.floor(this.lower[i] + 0.5);
            int up = (int) Math.floor(this.upper[i] + 0.5);
            //IloNumVar xi = this.cplex.intVar(low, up, varname);
            IloNumVar xi = this.cplex.intVar(low, Math.max(low, up), varname);
            this.x.add( xi);
        }
        for ( int i = this.numberOfIntegerVariables ; i< this.n; i++){
            String varname = "y" + i;
            //int low = (int) Math.floor(this.lower[i] + 0.5);
            //int up = (int) Math.floor(this.upper[i] + 0.5);
            //IloNumVar xi = this.cplex.intVar(low, up, varname);
            IloNumVar xi = this.cplex.numVar(lower[i], upper[i], varname);
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
            if(directions[i] > 0) {
                this.cplex.addGe(expr, b[i]);
            }else if( directions[i] < 0){
                this.cplex.addLe(expr, b[i]);
            }else{
                this.cplex.addEq(expr, b[i]);
            }
        }
    }




    private double[] fill(boolean up){
        double [] d = new double [this.n];
        for(int i = 0 ; i < n; i++){
            d[i] = (up? Double.MAX_VALUE : -(Double.MAX_VALUE -1));
        }
        return d;
    }


    private ArrayList<Double> getX() throws IloException{
        ArrayList<Double> x = new ArrayList<>();
        for(int i = 0; i < this.n; i++){
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

            if(directions[i] > 0) {
                if (val < b[i] - eps) return false;
            }else if( directions[i] < 0){
                if (val > b[i] + eps) return false;
            }else{
                if (val < b[i] - eps  || val > b[i] + eps  ) return false;
            }
        }
        return true;
    }



    public ArrayList<Object> solve() throws IloException{
        this.cplex.setOut(null);
        //this.cplex.setParam(IloCplex.DoubleParam.EpGap, 0d);
        //this.cplex.setParam(IloCplex.Param.MIP.Strategy.FPHeur, 1);
        this.cplex.setParam(IloCplex.Param.MIP.Limits.Solutions, 2);

        long start = System.currentTimeMillis();
        //System.out.println(this.model);
        this.cplex.solve();
        //
        ArrayList<Object> ret = new ArrayList<>();
        ArrayList<Double> x = this.getX();
            long end = System.currentTimeMillis();
            ret.add(x);
            ret.add(true);
            ret.add(0);
            ret.add(-start + end);
            ret.add(this.checkFeasibility(x));
        ret.add( FButils.objVal(x, this.c));
        return ret;
    }

}
