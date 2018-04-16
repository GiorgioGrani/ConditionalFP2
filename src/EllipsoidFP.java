import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.util.List;
import java.util.HashMap;
import java.util.Map;

import java.util.ArrayList;

public class EllipsoidFP  {
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
    private static final double safetythreshold = 5;

    private IloCplex cplex;
    private IloModel model;

    private IloAddable objective;
    private IloAddable distanceobj;
    private ArrayList< IloNumVar> x;
    private ArrayList<ArrayList<Double>> awarex;

    private double[] ellipsoid ;
    private ArrayList<Double> center;
    private double rayq;




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
        this.sign = Math.pow(-1d, this.numberOfIntegerVariables + 1);//*Math.pow(2,-this.numberOfIntegerVariables);


        double max = -(Double.MAX_VALUE - 1d);
        double min = Double.MAX_VALUE;
        for(int i = 0 ; i < this.numberOfIntegerVariables; i++) {
            if(c[i] > max) max = c[i];
            if(c[i] < min) min = c[i];
        }
        this.ellipsoid = new double[this.numberOfIntegerVariables];
        for(int i = 0 ; i < this.ellipsoid.length ; i++){
            this.ellipsoid[i] = (max+1 - c[i])/(max-min+1);
            System.out.println("elipsoid" +ellipsoid[i]);
        }


        ArrayList<Double> center = new ArrayList<>();
        double ray = 0;
        for(int i = 0 ; i < this.numberOfIntegerVariables; i++) {
            double val = this.upper[i] - this.lower[i];
            ray += Math.pow(val, 2);
            center.add((this.c[i]>0? this.upper[i] : this.lower[i]));
        }
        this.rayq = ray + this.safetythreshold;
        this.center = center;

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

    private void setDistanceObjective( ArrayList<Double> grad0) throws IloException{
        IloLinearNumExpr expr = this.cplex.linearNumExpr();
        double maxabs = 0d;
        for(Double d : grad0){
            if(Math.abs(d)>maxabs) maxabs = Math.abs(d);
        }
        ArrayList<Double> grad = new ArrayList<>();
        System.out.println(".|.|.|    "+maxabs);
        if(maxabs < 1 && maxabs > 0){
            for(Double d : grad0){
                grad.add(d/maxabs);
                //System.out.println("mod grad"+(d/maxabs));
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

    private ArrayList<Double> getDirection(ArrayList<Double> x, ArrayList<Double> xk){
        ArrayList<Double> dk = new ArrayList<>();
        int i = 0;
        for(Double xi : x){
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
            double val = xi + ak*dk.get(i);
            step.add(val);
            i++;
        }
        return step;
    }

    private ArrayList<Double> nextStepArmijo(ArrayList<Double> xk, double ak, ArrayList<Double> dk){
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

    private boolean checkFlatZone(ArrayList<Double> xk){
        int i = 0;
        for(Double xi : xk){
            if( Math.abs(xi - Math.floor(xi + 0.5)) >= (0.5 - eps)) {
                return true;
            }
            i++;
            if(i>=this.numberOfIntegerVariables) break;
        }
        return false;
    }

    private double Armijo(ArrayList<Double> xk, ArrayList<Double> dk, ArrayList<Double> grad){
        double ak = 1;
        double g = 0.001;
        double d = 0.9;
        double value = scalarProd(dk, grad);
        double norm = FButils.Chebynorm(dk,this.numberOfIntegerVariables);
        if(norm > 1) {
            ak = 1d/norm;
            ArrayList<Double> xk1 = nextStepArmijo(xk, ak, dk);

            if (f(xk) >= 1d+ellipsoidf(xk)[0] && f(xk1) <1d+ellipsoidf(xk)[0] && !checkFlatZone(xk1)) {
                System.out.println("_______________--------------------- TADAAAAAAAAN");
                return ak;
            }
        }
        System.out.println(" Armijo-Control "+checkFlatZone(nextStepArmijo(xk, ak, dk)));
        ArrayList<Double> xk1 = nextStepArmijo(xk, ak, dk);
        while( f(xk1) > f(xk) + g*ak*value || checkFlatZone(xk1)){
            System.out.println(norm+"Armi]]]]jo-Control   "+f(nextStep(xk,ak,dk))+"  ?<= " +(f(xk) + g*ak*value)+" ak "+ak);
            ak = ak*d;
            xk1 = nextStepArmijo(xk, ak, dk);
        }
        return ak;
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

    private ArrayList<Double> function(ArrayList<Double> x){
        double ret = 1;
        ArrayList<Double> res = new ArrayList<>();
        int i = 0;
        //System.out.println("RET "+checkFlatZone(x));
        for(Double xi : x){
            if(i >= this.numberOfIntegerVariables) break;

            ret = ret*(-1d-Math.cos(2*Math.PI*xi))/2d ;
            //if(i >= this.numberOfIntegerVariables)System.out.print("("+xi+" , "+ret+")");
            i++;
        }
        System.out.println();

        ret = this.sign*ret;
        double func = Math.exp(ret);
        res.add(func);
        res.add(ret);
        return res;
    }

    private double f(ArrayList<Double> x){
        double ff = this.filledfunction(x);
        if(ff > 0){

            return ff + ellipsoidf(x)[0];
        }
        //System.out.println("Funcfunc "+ this.function(x).get(0));
        return this.function(x).get(0) + ellipsoidf(x)[0];
    }

    private double[] ellipsoidf(ArrayList<Double> x){
        double ret = 0;
        int i = 0;
        for(Double xi : x){
            ret += Math.pow(xi - this.center.get(i) , 2 )*this.ellipsoid[i];
            //System.out.println("               ejejeje   ) "+this.center.get(i));

            i++;
            if(i >= this.numberOfIntegerVariables) break;
        }
        double [] res = new double [2];
        res[0] = - Math.log(rayq-ret);
        res[1] = rayq - ret;
        return res;
    }

    private double filledfunction(ArrayList<Double> x){
        if(this.awarex == null) return 0d;
        for(ArrayList<Double> xh : this.awarex){
            if(FButils.Chebynorm(xh,x, this.numberOfIntegerVariables) <= 0.5){
                //System.out.println("------------------------------------"+(2d - function(x).get(0)));
                return 1d;
            }
        }
        return 0d;
    }


    private ArrayList<Double> grad(ArrayList<Double> x){
        boolean filled = this.filledGrad(x);
        if (filled) return  ellipsoidg(x);
        return FButils.sum(this.functiongrad(x), ellipsoidg(x), this.numberOfIntegerVariables);
    }

    private ArrayList<Double> ellipsoidg(ArrayList<Double> x ){
        ArrayList<Double> ret = new ArrayList<>();
        double div = ellipsoidf(x)[1];
        int i = 0;
        for(Double xi : x){
            double num =  2d*ellipsoid[i]*(xi - this.center.get(i))/div;
            //System.out.println("               ;;;   ) "+num);
            ret.add(num);
            i++;
            if(i >= this.numberOfIntegerVariables) break;
        }
        return ret;
    }

    private ArrayList<Double> functiongrad(ArrayList<Double> x){
        ArrayList<Double> ret = new ArrayList<Double>();
        ArrayList<Double> singlevals = new ArrayList<>();
        int i = 0;
        for(Double xi : x){
            if(i >= this.numberOfIntegerVariables) break;
            singlevals.add((-1d-Math.cos(2*Math.PI*xi))/2d) ;
            i++;
        }
        i = 0;
        ArrayList<Double> res = function(x);
        double func = res.get(0);
        double px = res.get(1);
        double sign = 1d;


        for(Double xi : x){
            if(i >= this.numberOfIntegerVariables) break;
            double val = sign*Math.PI*Math.sin(2*Math.PI*xi);
            //System.out.print(px+"__"+val);
            if(Math.abs(singlevals.get(i)) >= 1e-20){
                val = val*px/singlevals.get(i);

                val = val*func;
            }else{
                val = 0d;
            }
            //System.out.println("++"+val);
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
            this.awarex = new ArrayList<>();
            this.awarex.add(xh);
            return;
        }

        this.awarex.add(xh);
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
            ret.add( FButils.objVal(xk, this.c));
            return ret;
        }

        boolean print = false;


        int count = 0;
        this.model.remove(this.objective);
        while(!stop && maxiter <= 10) {
            System.out.println("ITER>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "+maxiter);
            maxiter++;
            //print(xk);
            //System.out.println("BOIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
            //print(grad);
            this.setDistanceObjective( grad);
            System.out.println("argh");
            this.cplex.solve();
//            if((maxiter-1) % 100 == 0)
            System.out.println(maxiter+") "+FButils.L2Norm(grad));
//            if(FButils.L2Norm(grad) <= 1e-2){
//                for(int i = 0; i < this.numberOfIntegerVariables; i++){
//                    if(Math.random() > 0.8)
//                    this.center.set(i, this.center.get(i)+(2*Math.random()-1d));
//                }
//            }
            //double val = scalarProd(grad,x) - scalarProd(grad, xk);

//if(Math.abs(val) <= 1e-10 || f(xk) < 1d){
            System.out.println(f(xk)+"  "+ellipsoidf(xk)[0]+"  "+ (f(xk)<1d));
            if( f(xk) < 1d + ellipsoidf(xk)[0]){
                count++;
                System.out.println("fillato");
                ArrayList<Double> xh = this.getXTilde(xk);


                this.model.remove(this.distanceobj);
                this.model.add(this.objective);
                for(int i = 0; i < this.numberOfIntegerVariables; i++){
                    double value = xh.get(i);
                    this.x.get(i).setLB(value);
                    this.x.get(i).setUB(value);
                }



                feas = this.cplex.solve();
                if(!feas){
                    //System.out.println("   ahio");
                    if(print) System.out.println("                                    filling"+FButils.L2Norm(grad));
                    if(print) print(xh);
                    if(print) System.out.println("                    KK");
                    this.fill(xh);
                    if(print) print(grad);
                    if(print) System.out.println("                   KKGRADBEFOREFILLING<<");
                    grad = this.grad(xk);
                    if(print) print(grad);
                    if(print) System.out.println("                   KKGRADAFTER");

                }else{

                    for(int i = this.numberOfIntegerVariables; i < this.n; i++){
                        xh.add(this.cplex.getValue(this.x.get(i)));
                    }
                    xk = xh;
                    stop = true;
                }
                for(int i = 0; i < this.numberOfIntegerVariables; i++){

                    this.x.get(i).setLB(this.lower[i]);
                    this.x.get(i).setUB(this.upper[i]);
                }
                this.model.remove(this.objective);
            }else{
                if(print) System.out.println("obj:" + this.cplex.getModel());
                if(print) System.out.println("------------------------xk");
                if(print) print(xk);
                if(print) System.out.println("------------------------");

                x = this.getX();
                if(print) System.out.println("------------------------xx");
                if(print) print(x);
                if(print) System.out.println("------------------------(( "+stopCondition(x));

                ArrayList<Double> dk = this.getDirection(x, xk);
                if(print) System.out.println("------------------------DDDD");
                if(print) print(dk);
                if(print) System.out.println("------------------------");

                double ak = this.Armijo(xk, dk, grad);
                xk = this.nextStep(xk, ak, dk);
                if(print) System.out.println("------------------------xk+1");
                if(print) print(xk);
                if(print) System.out.println("------------------------<<");
                System.out.println("is here");
                stop = this.stopCondition(xk);
                System.out.println("                "+checkFeasibility(xk));
                if(!stop) {

                    grad = this.grad(xk);
                }
                //             print(grad);
                //           System.out.println("                    KKGRADBEFORE  ak "+ak);
                //         System.out.println(maxiter+">> "+(-val)+" <= "+1e-10+"     "+(Math.abs(val) <= 1e-10)+"  f "+f(xk));
            }
            //System.out.println(maxiter+") "+val +"    "+FButils.integralityGap(IntegralityGapTypes.L2Norm,this.numberOfIntegerVariables,xk));
            this.model.remove(this.distanceobj);
        }
        long end = System.currentTimeMillis();
        print(xk);
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
        int i = 0;
        for(Double xi : xk){
            System.out.println("                         "+xi+" "+this.center.get(i));
            i++;
            if(i>= this.numberOfIntegerVariables) break;
        }
    }

}

