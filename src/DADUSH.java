import ilog.concert.*;
import ilog.cplex.IloCplex;
import java.util.Collections;
import java.util.Random;

import java.util.ArrayList;

public class DADUSH  {
    private double[] c;
    private double[][] A;
    private double [][] V;
    private double [] b;
    private double [] lower;
    private double [] upper;
    private int [] directions;
    private int numberOfIntegerVariables;
    private int n ;
    private double sign;
    private static final double eps = 1e-6;
private Random rand;

    private IloCplex cplex;
    private IloModel model;

    private IloAddable objective;
    private IloAddable distanceobj;
    private ArrayList<IloNumVar> x;

    public void set(double[] c, double[][] A , double [] b, double [] lower , double [] upper, int numberOfIntegerVariables, int [] directions) throws IloException {
        this.rand = new Random(3);
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

        this.V = new double [A[0].length][A.length];
        double [] norms = new double [A.length];
        for( int i = 0 ; i < this.A.length ; i++){
            int j = 0;
            double norm = 0;
            IloLinearNumExpr expr = this.cplex.linearNumExpr();
            for( IloNumVar xi : x){
                if(Math.abs(A[i][j]) > 0) {
                    expr.addTerm(A[i][j], xi);
                    norm += A[i][j]*A[i][j];
                }
                j++;
            }
            if(directions[i] > 0) {
                this.cplex.addGe(expr, b[i]);
            }else if( directions[i] < 0){
                this.cplex.addLe(expr, b[i]);
            }else{
                this.cplex.addEq(expr, b[i]);
            }
            norms[i] = Math.sqrt(norm/(b[i]*b[i]));
        }

        double nmax = 0;
        for(double d : norms){
            if( nmax < d) nmax = d;
        }

        for(int i = 0 ; i < A.length; i++){
            for(int j = 0; j < A[0].length; j++){
                if(directions[i] > 0) {
                    V[j][i] = -A[i][j]/(b[i]*nmax);
                }else if( directions[i] < 0){
                    V[j][i] = A[i][j]/(b[i]*nmax);
                }else{
                    V[j][i] = A[i][j]/(b[i]*nmax);
                }
            }
        }

    }


    private double[] GSwalk(ArrayList<Integer> A){
        Collections.sort(A);
        //double [] ret = new double [this.numberOfIntegerVariables];
        double [][] B = new double [V.length][V[0].length];
        ArrayList<Integer> subA = new ArrayList<>();
        for(Integer i : A){
            B[i] = V[i];
            for(Integer j : subA){
               double mu = 0;
               mu = scalarProd(B[j], B[i]);
               double normb = FButils.L2Norm(B[j]);
               B[i] = update(B[i], mu/normb, B[j]);
            }
            subA.add(i);
        }

        return B[B.length-1];
    }

    private double [] update(double [] vb, double mu, double [] vj){
        double [] ret = vb;
        for(int i = 0  ; i < vb.length; i++){
            ret[i] = ret[i] - mu*vj[i];
        }
        return ret;
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
    private double scalarProd(double [] a , double [] b){
        double ret = 0;
        for(int i = 0; i < a.length; i++){
           ret += a[i]*b[i];
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



    private ArrayList<Integer> setA( ArrayList<Double> xk){
        ArrayList<Integer> ret = new ArrayList<>();
        int i = 0;
        //System.out.println();
        for(Double d : xk){
            //System.out.println("--> " + d);
            if(Math.abs(Math.round(d)-d) > 0){

                ret.add(i);
            }
            i++;
        }
        return ret;
    }

    private double [] deltaVal( ArrayList<Double> xk,double []  weights, ArrayList<Integer> A){
        double dm = Double.MAX_VALUE -1;
        double dp =  - (Double.MAX_VALUE -1 );

        int i = 0;
        for(Integer ind : A){
            double x = xk.get(ind);
            double val = weights[i];
            i++;
if(Math.abs(val) > 0) {
    double d1 = -2d * x / val;
    double d2 = (2d - 2d * x) / val;

    if (dm > d1 || dm > d2) {
        dm = Math.min(d1, d2);
    }
    if (dp < d1 || dp < d2) {
        dp = Math.max(d1, d2);
    }
}
        }

        double [] ret = new double [2];
        ret[0] = dm;
        ret[1] = dp;
        //System.out.println(dp+" "+dm);
        return ret;
    }

    private double [] diff(double [] a , double [] b){
        double [] ret = new double [a.length];
        for(int i = 0; i < a.length; i++){
            ret[i] = a[i] - b[i];
        }
        return ret;
    }

    private double [] solveSystem(double [] vnt,ArrayList<Integer> A) throws IloException{
        int n = A.size();
        //System.out.println(n);
        double [] rif = diff(vnt, V[A.get( n-1)]);


        IloCplex sub = new IloCplex();
        sub.setOut(null);
        //sub.setParam(IloCplex.Param.Preprocessing.Presolve, false);
        ArrayList<IloNumVar> vars = new ArrayList<>();
        for(int i = 0; i < n-1; i++){
            vars.add( sub.numVar(-(Double.MAX_VALUE-1), Double.MAX_VALUE));
        }

        for(int i = 0; i < vnt.length; i++){
            IloLinearNumExpr expr = sub.linearNumExpr();
            int ind = 0;
            for(Integer j : A){
                if(j.intValue() != A.get(n-1).intValue()){
                    //System.out.println(n+" "+ind);
                    expr.addTerm(V[j][i], vars.get(ind));
                }
                ind ++;
            }
            sub.addGe(expr, rif[i]-this.eps);sub.addLe(expr, rif[i]+this.eps);
        }

//        int ind = 0;
//        for(Integer i : A){
//            if( i.intValue() != A.get(n-1).intValue()) {
//                double val = 0;
//                for (int j = 0; j < V[0].length; j++) {
//                    val += V[i][j];
//                }
//                System.out.println("val "+ val+"  i " +i);
//                if (val <= 0) {
//                    System.out.println("Verificato " + ind+" "+i);
//                    IloLinearNumExpr expr = sub.linearNumExpr();
//                    expr.addTerm(vars.get(ind),1d);
//                    sub.addEq(expr, 0d);
//                }
//            }
//            ind ++;
//        }

        IloLinearNumExpr obj = sub.linearNumExpr();
        for(IloNumVar v : vars){
            obj.addTerm(0d, v);
        }
        sub.addMinimize(obj);
        sub.solve();
        System.out.println(sub.getModel());

        double []  ret = new double [n];
        int i = 0;
        for(IloNumVar v : vars){
            ret[i] = sub.getValue(v);
            i++;
        }
        ret[n-1] = 1d;
        return ret;


    }

    private ArrayList<Double> update(ArrayList<Double> xk, double [] delta, double [] weights, ArrayList<Integer> A){
        double dm = delta[0];
        double dp = delta[1];
        double probm = dp/(dp + Math.abs(dm));

        int i = 0;
        double set = 0;
        for(Integer ind : A){
            if( rand.nextDouble() < probm){
                set = dm;
            }else{
                set = dp;
            }
            xk.set(ind, xk.get(ind) + set+weights[i] );
            i++;
        }
        return xk;
    }

    private ArrayList<Double> normalize( ArrayList<Double> xk){
        double max = 0;
        double min = 0;
        for(Double xi : xk){
            if(xi > max) max = xi;
            if(xi < min) min = xi;
        }
        ArrayList<Double> newx = new ArrayList<>();

        double c = max - min;
        double rand = this.rand.nextDouble()*c+this.eps;
        c = c+ rand;
        //System.out.println(c);
        for(Double xi : xk){
            newx.add((xi - min)/c);
        }
        return newx;
    }

    public ArrayList<Object> solve() throws IloException{
        int maxiter = 0;
        this.cplex.setOut(null);
        this.cplex.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Barrier);
        //this.cplex.setParam(IloCplex.Param.Barrier.Limits.Iteration, 10);
        long start = System.currentTimeMillis();
        boolean solve = this.cplex.solve();
        //System.out.println(this.model);
        ArrayList<Object> ret = new ArrayList<>();
        ArrayList<Double> x = this.getX();
        ArrayList<Double> xk = x;

        boolean stop = this.stopCondition(xk);
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

        ArrayList<Integer> A = setA(xk);

        while(!stop && maxiter < 100){
            System.out.println(maxiter);
            maxiter++;
            double [] vnt = GSwalk(A);
            double [] weights;
            weights = solveSystem(vnt,A);
            double [] delta ;
            delta = deltaVal(xk, weights, A);
            xk = update(xk, delta, weights, A);
            xk = normalize(xk);
            xh = getXTilde(xk);
            if(checkFeasibility(xh)){
                xk = xh;
                stop = true;
                //System.out.println("cane");
            }

            A = setA(xk);

        }

//        while(!stop && maxiter <= 50) {
//            maxiter++;
//            this.setDistanceObjective( G);
//            this.cplex.solve();
//
//            //System.out.println("Here I am "+ count+" "+val+" "+(Math.abs(val) <= 1e-10)+" "+f(xk));
//            xh = this.getXTilde(xk);
//
//
//            this.model.remove(this.distanceobj);
//            this.model.add(this.objective);
//            for(int i = 0; i < this.numberOfIntegerVariables; i++){
//                double value = xh.get(i);
//                this.x.get(i).setLB(value);
//                this.x.get(i).setUB(value);
//            }
//            feas = this.cplex.solve();
//            if(feas){
//                for(int i = this.numberOfIntegerVariables; i < this.n; i++){
//                    xh.add(this.cplex.getValue(this.x.get(i)));
//                }
//                xk = xh;
//                stop = true;
//                break;
//            }
//            for(int i = 0; i < this.numberOfIntegerVariables; i++){
//                this.x.get(i).setLB(this.lower[i]);
//                this.x.get(i).setUB(this.upper[i]);
//            }
//            this.model.remove(this.objective);
//        }





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
            Double rd = FButils.randomRound(xi, this.c[i], this.rand);
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

