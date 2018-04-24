import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.util.ArrayList;

public class BasicVersion  {
    private double[] c;
    private double[][] A;
    private double [] b;
    private double [] lower;
    private double [] upper;
    private int [] directions;
    private int numberOfIntegerVariables;
    private int n ;
    private double sign;
    private static final double eps = 1e-6;

    private IloCplex cplex;
    private IloModel model;

    private IloAddable objective;
    private IloAddable distanceobj;
    private ArrayList<IloNumVar> x;

    public void set(double[] c, double[][] A , double [] b, double [] lower , double [] upper, int numberOfIntegerVariables, int [] directions) throws IloException {
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

        for(int i = 0; i < this.numberOfIntegerVariables; i++){
            this.lower[i] = 0d;
            this.upper[i] = 1d;

            //todo hai supposto solo binarie come parte intera
        }

        this.directions = directions;
        this.cplex = new IloCplex();
        this.model = this.cplex.getModel();
        this.sign = Math.pow(-1d, this.numberOfIntegerVariables + 1);//*Math.pow(2,-this.numberOfIntegerVariables);

        this.createVariables();
        this.setObjective();
        this.setConstraints();
    }

    private double[] fill(boolean up){
        double [] d = new double [this.n];
        for(int i = 0 ; i < n; i++){
            d[i] = (up? Double.MAX_VALUE : -(Double.MAX_VALUE -1));
        }
        return d;
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
            if(directions[i] > 0) {
                this.cplex.addGe(expr, b[i]);
            }else if( directions[i] < 0){
                this.cplex.addLe(expr, b[i]);
            }else{
                this.cplex.addEq(expr, b[i]);
            }
        }
    }







    private ArrayList<Double> getX() throws IloException{
        ArrayList<Double> x = new ArrayList<>();
        for(int i = 0; i < this.n; i++){
            Double d = this.cplex.getValue(this.x.get(i));
            x.add(d);
        }
        return x;
    }

    private double g(double d){
        if (Math.abs(d) <= this.eps) return 1d;
        if (Math.abs(1-d) <= this.eps) return -1d;
        return 0d;
    }

    private ArrayList<Double> G(ArrayList<Double> x0){
        ArrayList<Double> ret = new ArrayList<>();
        int i = 0;
        for(Double d: x0){
            ret.add(g(d));
            i++;
            if(i >= this.numberOfIntegerVariables) return ret;
        }
        return ret;
    }

    private void setDistanceObjective( ArrayList<Double> Gx0) throws IloException{
        IloLinearNumExpr expr = this.cplex.linearNumExpr();

        int i = 0;
        for( Double d : Gx0){
            //double d = grad.get(i);
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


    private double scalarProd(ArrayList<Double> a, ArrayList<Double> b){
        double ret = 0;
        int i = 0;
        if(a.size() > b.size()){
            for (Double ai : b) {
                ret = ret + ai * a.get(i);
                i++;
            }
        }else {
            for (Double ai : a) {
                ret = ret + ai * b.get(i);
                i++;
            }
        }
        return ret;
    }
    private double scalarProdInt(ArrayList<Double> a, ArrayList<Double> b){
        double ret = 0;
        int i = 0;
        if(a.size() > b.size()){
            for (Double ai : b) {
                ret = ret + ai * a.get(i);
                i++;
                if(i >= this.numberOfIntegerVariables) return ret;
            }
        }else {
            for (Double ai : a) {
                ret = ret + ai * b.get(i);
                i++;
                if(i >= this.numberOfIntegerVariables) return ret;
            }
        }
        return ret;
    }




    public ArrayList<Object> solve() throws IloException{
        int maxiter = 0;
        this.cplex.setOut(null);
        //this.cplex.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Barrier);

        long start = System.currentTimeMillis();
        this.cplex.solve();
        //System.out.println(this.model);
        ArrayList<Object> ret = new ArrayList<>();
        ArrayList<Double> x = this.getX();
        ArrayList<Double> xk = x;

        boolean stop = this.stopCondition(xk);
        boolean feas = false;
        double lpobj = this.cplex.getObjValue();
        if(stop){
            long end = System.currentTimeMillis();
            ret.add(x);
            ret.add(stop);
            ret.add(0);
            ret.add(-start + end);
            ret.add(this.checkFeasibility(x));
            ret.add( FButils.objVal(xk, this.c));
            return ret;
        }

        boolean print = false;

        ArrayList<Double> xh = this.getXTilde(xk);
        ArrayList<Double> G = G(xh);
        this.model.remove(this.objective);
        while(!stop && maxiter <= 50) {
            maxiter++;
            this.setDistanceObjective( G);
            this.cplex.solve();

            //System.out.println("Here I am "+ count+" "+val+" "+(Math.abs(val) <= 1e-10)+" "+f(xk));
            xh = this.getXTilde(xk);


            this.model.remove(this.distanceobj);
            this.model.add(this.objective);
            for(int i = 0; i < this.numberOfIntegerVariables; i++){
                double value = xh.get(i);
                this.x.get(i).setLB(value);
                this.x.get(i).setUB(value);
            }
            feas = this.cplex.solve();
            if(feas){
                for(int i = this.numberOfIntegerVariables; i < this.n; i++){
                    xh.add(this.cplex.getValue(this.x.get(i)));
                }
                xk = xh;
                stop = true;
                break;
            }
            for(int i = 0; i < this.numberOfIntegerVariables; i++){
                this.x.get(i).setLB(this.lower[i]);
                this.x.get(i).setUB(this.upper[i]);
            }
            this.model.remove(this.objective);
        }
        long end = System.currentTimeMillis();
        //print(xk);
        ret.add(xk);
        ret.add(stop);
        ret.add(maxiter);
        ret.add(-(start - end));
        ret.add(this.checkFeasibility(xk));
        ret.add( FButils.objVal(xk, this.c));
        return ret;
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

    private void print(ArrayList<Double> xk){
        for(Double xi : xk){
            System.out.println("                         "+xi);
        }
    }

}

