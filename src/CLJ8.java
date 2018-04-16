import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.util.List;
import java.util.HashMap;
import java.util.Map;

import java.util.ArrayList;

public class CLJ8  {
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
    private ArrayList<ArrayList<Double>> awarex;

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
    private ArrayList<Double> getX() throws IloException{
        ArrayList<Double> x = new ArrayList<>();
        for(int i = 0; i < this.numberOfIntegerVariables; i++){
            Double d = this.cplex.getValue(this.x.get(i));
            x.add(d);
        }
        return x;
    }

    private void setDistanceObjective( ArrayList<Double> grad0) throws IloException{
        IloLinearNumExpr expr = this.cplex.linearNumExpr();
        double maxabs = 0d;
        for(Double d : grad0){
            if(Math.abs(d)>maxabs) maxabs = Math.abs(d);
        }
        ArrayList<Double> grad = new ArrayList<>();
        if(maxabs < 1 && maxabs > 0){
            for(Double d : grad0){
                grad.add(d/maxabs);
                // System.out.println("mod grad"+(d/maxabs));
            }
        }else{
            grad = grad0;
        }
        int i = 0;
        for( Double d : grad){
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

    private double Armijo(ArrayList<Double> xk, ArrayList<Double> dk, ArrayList<Double> grad){
        double ak = 1;
        double g = 0.01;
        double d = 0.5;

        double value = scalarProd(dk, grad);
        while( f(nextStep(xk,ak,dk)) > f(xk) + g*ak*value){
            // System.out.println("Armijo-Control   "+f(nextStep(xk,ak,dk))+"  ?<= " +(f(xk) + g*ak*value+" ak "+ak));
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

    private ArrayList<Double> function(ArrayList<Double> x){
        double ret = 1;
        ArrayList<Double> res = new ArrayList<>();
        for(Double xi : x){
            ret = ret*(-1d-Math.cos(2*Math.PI*xi)) ;
        }

        ret = Math.pow(-1d, this.numberOfIntegerVariables + 1)*Math.pow(2,-this.numberOfIntegerVariables)*ret;
        double func = (ret);
        res.add(func);
        res.add(ret);
        return res;
    }

    private double f(ArrayList<Double> x){
        double ff = this.filledfunction(x);
        if(ff > 0){
            //System.out.println("Funcff "+ff);
            return ff;
        }
        //System.out.println("Funcfunc "+ this.function(x).get(0));
        return this.function(x).get(0);
    }

    private double filledfunction(ArrayList<Double> x){
        if(this.awarex == null) return 0d;
        for(ArrayList<Double> xh : this.awarex){
            if(FButils.Chebynorm(xh,x, this.numberOfIntegerVariables) <= 0.5){
                //System.out.println("------------------------------------"+(2d - function(x).get(0)));
                return  - function(x).get(0);
            }
        }
        return 0d;
    }


    private ArrayList<Double> grad(ArrayList<Double> x){
        boolean filled = this.filledGrad(x);

        return this.functiongrad(x, filled);
    }

    private ArrayList<Double> functiongrad(ArrayList<Double> x, boolean filled){
        ArrayList<Double> ret = new ArrayList<Double>();
        ArrayList<Double> singlevals = new ArrayList<>();
        for(Double xi : x){
            singlevals.add(-1d-Math.cos(2*Math.PI*xi)) ;
        }
        int i = 0;
        ArrayList<Double> res = function(x);
        double func = res.get(0);
        double px = res.get(1);
        double sign = 1d;
        if(filled) sign = -1d;

        for(Double xi : x){
            double val = sign*2d*Math.PI*Math.sin(2*Math.PI*xi);
            val = val*px/singlevals.get(0);
            //val = val*func;
            ret.add(val);
            i++;
        }
        return ret;
    }

    private boolean filledGrad(ArrayList<Double> x){
        if(this.awarex == null){
            return false;
        }
        for(ArrayList<Double> xh : this.awarex){
            if(FButils.Chebynorm(xh,x, this.numberOfIntegerVariables) <= 0.5){
                return true;
            }
        }
        return false;
    }


    private void fill(ArrayList<Double> xh){
        if(this.awarex == null){
            this.awarex = new ArrayList<ArrayList<Double>>();
            this.awarex.add(xh);
            return;
        }

        this.awarex.add(xh);
    }

    public ArrayList<Object> solve() throws IloException{
        int maxiter = 0;
        this.cplex.setOut(null);

        long start = System.currentTimeMillis();
        this.cplex.solve();
        ArrayList<Object> ret = new ArrayList<>();
        ArrayList<Double> x = this.getX();
        ArrayList<Double> xk = x;
        ArrayList<Double> grad = this.grad(xk);
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
            ret.add( this.checkGap(x,lpobj));
            return ret;
        }



        int count = 0;
        this.model.remove(this.objective);
        while(!stop && maxiter <= 5) {
//System.out.println("ITER>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "+maxiter);
            maxiter++;
            this.setDistanceObjective( grad);
            this.cplex.solve();
            System.out.println(maxiter+") "+FButils.L2Norm(grad));

            double val = scalarProd(grad,x) - scalarProd(grad, xk);

//if(Math.abs(val) <= 1e-10 || f(xk) < 1d){

            if( f(xk) < 0d){
                count++;
                //System.out.println("Here I am "+ count+" "+val+" "+(Math.abs(val) <= 1e-10)+" "+f(xk));
                ArrayList<Double> xh = this.getXTilde(xk);

                feas = this.checkFeasibility(xh);
                if(!feas){
                    //  System.out.println("                                    filling"+FButils.L2Norm(grad));
                    //print(xh);
                    //System.out.println("                    KK");
                    this.fill(xh);
                    //print(grad);
                    //System.out.println("                   KKGRADBEFOREFILLING<<");
                    grad = this.grad(xk);
                    //print(grad);
                    //System.out.println("                   KKGRADAFTER");
                }else{
                    xk = xh;
                    stop = true;
                }

            }else{
                //System.out.println("obj:" + this.cplex.getModel());
                //System.out.println("------------------------xk");
                //print(xk);
                //System.out.println("------------------------");

                x = this.getX();
                //System.out.println("------------------------xx");
                //print(x); System.out.println("------------------------(( "+stopCondition(x));

                ArrayList<Double> dk = this.getDirection(x, xk);
                //System.out.println("------------------------DDDD");
                //print(dk);
                //System.out.println("------------------------");

                double ak = this.Armijo(xk, dk, grad);
                xk = this.nextStep(xk, ak, dk);
                //System.out.println("------------------------xk+1");
//print(xk); System.out.println("------------------------<<");

                stop = this.stopCondition(xk);
                if(!stop) grad = this.grad(xk);
                //             print(grad);
                //           System.out.println("                    KKGRADBEFORE  ak "+ak);
                //         System.out.println(maxiter+">> "+(-val)+" <= "+1e-10+"     "+(Math.abs(val) <= 1e-10)+"  f "+f(xk));
            }
            //System.out.println(maxiter+") "+val +"    "+FButils.integralityGap(IntegralityGapTypes.L2Norm,this.numberOfIntegerVariables,xk));
            this.model.remove(this.distanceobj);
        }
        long end = System.currentTimeMillis();
        ret.add(xk);
        ret.add(stop);
        ret.add(maxiter-1);
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
            System.out.println("                       - "+xi);
        }
    }

}

