import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.util.ArrayList;
import java.util.Collections;

public class CDFP  {
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
    private ArrayList<ArrayList<Double>> points = new ArrayList<>();

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

    private void setDistanceObjective( double[] Gx0) throws IloException{
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

    private ArrayList<Ord> giveOrder() {

        ArrayList<Ord> ord = new ArrayList<>();
        double norm = 0d;
        for (int j = 0; j < this.c.length ; j++) {
            double xi = c[j];
            norm += xi * xi;
        }
        norm = Math.sqrt(norm);


        for (int i = 0; i < this.A.length; i++) {
            double val = 0;
            double innernorm = 0d;
            for (int j = 0; j < this.A[0].length; j++) {
                double xi = c[j];
                val += xi * A[i][j];
                innernorm += A[i][j]*A[i][j];
            }

            innernorm = Math.sqrt(innernorm);

            if (directions[i] < 0) {
                val = -val;
            }

            if(Math.abs(directions[i]) > 0 ){
                ord.add(new Ord( i , val/(norm*innernorm)));
            }
        }

        Collections.sort(ord);
        return ord;
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

    private void reinforcementCut(ArrayList<Double> G, ArrayList<Double> xh, double reinforce) throws IloException{
        IloLinearNumExpr expr = this.cplex.linearNumExpr();
        //System.out.println();
        for(int i = 0; i < this.numberOfIntegerVariables; i++){
            expr.addTerm(G.get(i), this.x.get(i));
//            if(G.get(i) < 0d) {
//                System.out.print("( " + xh.get(i) + "  " + G.get(i) + " ) ");
//            }else{
//                System.out.print("                 ");
//            }
        }
        //System.out.println();
        //System.out.println(
        this.cplex.addGe(expr, reinforce);
        //);
    }


    private ArrayList<Object> checkProximality(ArrayList<Double> xh){
        return new ArrayList<Object>();
        //todo implementare
    }


    private void combinedReinforcementCut(ArrayList<Double> xh, double reinforce) throws IloException{
        ArrayList<Object> check = checkProximality(xh);
        boolean find = (boolean) check.get(0);
        if(!find) return;
        ArrayList<Double> xhm1 = (ArrayList<Double>) check.get(1);
        IloLinearNumExpr expr = this.cplex.linearNumExpr();
        System.out.println();
        for(int i = 0; i < this.numberOfIntegerVariables; i++){
            double d = g(xh.get(i))+g(xhm1.get(i));
            expr.addTerm(d, this.x.get(i));
        }
        //System.out.println();
        System.out.println(this.cplex.addGe(expr, reinforce));

    }

    private ArrayList<Double> localSearch(ArrayList<Double> xh){
        ArrayList<Integer> indeces = new ArrayList<>();
        int i = 0;
        for(double xi : xh){
            if( this.c[i] > 0 && Math.abs(xi) <= this.eps){
                indeces.add(i);
            }else if( this.c[i] < 0 && Math.abs(xi) >= 1d-this.eps){
                indeces.add(i);
            }
            i++;
            if (i >= this.numberOfIntegerVariables) break;
        }
        if(indeces.get(0) == null) return xh;
        int minimumind = indeces.get(0);
        for(int ind : indeces){
            if(Math.abs(this.c[minimumind]) < Math.abs(this.c[ind])) minimumind = ind;
        }

        ArrayList<Double> ret = new ArrayList<>();
        int j = 0;
        for(Double xi : xh){
            double add = xi;
            if(j == minimumind){
                if(Math.abs(xi) >= 1d-this.eps) ret.add(0d);
                else ret.add(1d);
            }else{
                ret.add(xi);
            }
            j++;
            if(j >= this.numberOfIntegerVariables) break;
        }
        return ret;

    }


    public ArrayList<Object> solve() throws IloException {
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
        if (stop) {
            long end = System.currentTimeMillis();
            ret.add(x);
            ret.add(stop);
            ret.add(0);
            ret.add(-start + end);
            ret.add(this.checkFeasibility(x));
            ret.add(FButils.objVal(xk, this.c));
            return ret;
        }

        boolean print = false;

        ArrayList<Double> xh = this.getXTilde(xk);
        ArrayList<Double> xhm1 = new ArrayList<>();
        ArrayList<Double> G = G(xh);
        double valueobj = 1+scalarProdInt(G, xh);;
        this.model.remove(this.objective);
        int count = 0;
        int ind = 0;
        ArrayList<Ord> ord = this.giveOrder();
        while (!stop && ind < ord.size()) {
            maxiter++;
           // System.out.println("ci sono "+maxiter+"  "+ind);
//            this.reinforcementCut(G, xh, valueobj);
            this.setDistanceObjective(A[ord.get(ind).row]);
            ind ++;
            this.cplex.solve();


            xk = this.getX();


            //System.out.println("Here I am "+ count+" "+val+" "+(Math.abs(val) <= 1e-10)+" "+f(xk));
            xhm1 = xh;
            xh = this.getXTilde(xk);
            G = G(xh);
//            boolean proximity = FButils.L2norm(xhm1, xh, numberOfIntegerVariables) <= eps;
//            if (!proximity){
//                count = 0;
//                valueobj = 1d + scalarProdInt(G, xh);
//            } else if(proximity && count > 0){
//                count++;
//                valueobj ++;
//            }else if (proximity  && count ==0) {
//                double vobj = Math.ceil(this.cplex.getObjValue());
//                if(Math.abs(vobj - valueobj) <= this.eps){
//                    valueobj ++;
//                }else{
//                    valueobj = vobj;
//                }
//                count = 1;
//            }

            boolean equals = FButils.L2norm(xhm1, xh, this.numberOfIntegerVariables) <= eps;


//todo se sono uguali non ha senso riproposse cplex
            if (this.n - this.numberOfIntegerVariables > 0) {
                this.model.remove(this.distanceobj);
                this.model.add(this.objective);
                for (int i = 0; i < this.numberOfIntegerVariables; i++) {
                    double value = xh.get(i);
                    this.x.get(i).setLB(value);
                    this.x.get(i).setUB(value);
                }
                feas = this.cplex.solve();
                //System.out.println("feas  " + feas);
                if (feas) {
                    for (int i = this.numberOfIntegerVariables; i < this.n; i++) {
                        xh.add(this.cplex.getValue(this.x.get(i)));
                    }
                    xk = xh;
                    stop = true;
                    break;
                }


                for (int i = 0; i < this.numberOfIntegerVariables; i++) {
                    this.x.get(i).setLB(this.lower[i]);
                    this.x.get(i).setUB(this.upper[i]);
                }
                this.model.remove(this.objective);
            } else if (checkFeasibility(xh)) {
                xk = xh;
                stop = true;
                break;
            } else {
                this.model.remove(this.distanceobj);
            }


        }
        ind = 0;

        while (!stop && ind < this.numberOfIntegerVariables) {
            maxiter++;
//            this.reinforcementCut(G, xh, valueobj);
            this.setDistanceObjective(singleval(ind));
            ind ++;
            this.cplex.solve();


            xk = this.getX();


            //System.out.println("Here I am "+ count+" "+val+" "+(Math.abs(val) <= 1e-10)+" "+f(xk));
            xhm1 = xh;
            xh = this.getXTilde(xk);
            G = G(xh);
//            boolean proximity = FButils.L2norm(xhm1, xh, numberOfIntegerVariables) <= eps;
//            if (!proximity){
//                count = 0;
//                valueobj = 1d + scalarProdInt(G, xh);
//            } else if(proximity && count > 0){
//                count++;
//                valueobj ++;
//            }else if (proximity  && count ==0) {
//                double vobj = Math.ceil(this.cplex.getObjValue());
//                if(Math.abs(vobj - valueobj) <= this.eps){
//                    valueobj ++;
//                }else{
//                    valueobj = vobj;
//                }
//                count = 1;
//            }

            boolean equals = FButils.L2norm(xhm1, xh, this.numberOfIntegerVariables) <= eps;


//todo se sono uguali non ha senso riproposse cplex
            if (this.n - this.numberOfIntegerVariables > 0) {
                this.model.remove(this.distanceobj);
                this.model.add(this.objective);
                for (int i = 0; i < this.numberOfIntegerVariables; i++) {
                    double value = xh.get(i);
                    this.x.get(i).setLB(value);
                    this.x.get(i).setUB(value);
                }
                feas = this.cplex.solve();
                //System.out.println("feas  " + feas);
                if (feas) {
                    for (int i = this.numberOfIntegerVariables; i < this.n; i++) {
                        xh.add(this.cplex.getValue(this.x.get(i)));
                    }
                    xk = xh;
                    stop = true;
                    break;
                }


                for (int i = 0; i < this.numberOfIntegerVariables; i++) {
                    this.x.get(i).setLB(this.lower[i]);
                    this.x.get(i).setUB(this.upper[i]);
                }
                this.model.remove(this.objective);
            } else if (checkFeasibility(xh)) {
                xk = xh;
                stop = true;
                break;
            } else {
                this.model.remove(this.distanceobj);
            }


        }
        long end = System.currentTimeMillis();
        //print(x);
        ret.add(xk);
        ret.add(stop);
        ret.add(maxiter);
        ret.add(-(start - end));
        ret.add(this.checkFeasibility(xk));
        ret.add(FButils.objVal(xk, this.c));
        return ret;
    }
private double [] singleval(int ind){
        double [] ret = new double [this.numberOfIntegerVariables];
        ret[ind] = 1d;
        return ret;
}


    private ArrayList<Double> getXTilde(ArrayList<Double> x) throws IloException{
        ArrayList<Double> ret = new ArrayList<>();
        boolean stop = true;
        for(int i = 0; i < this.numberOfIntegerVariables; i++){
            Double xi = x.get(i);
            Double rd = FButils.round(xi);
            // Double rd = FButils.randomRound(xi, this.c[i]);

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

    private class Ord implements Comparable{
       public int row;
       public double value;

       public Ord(int i , double val){ row = i; value = val;}

       @Override
       public int compareTo(Object ob){
           Ord o = (Ord) ob;
           if(Math.abs(o.value) < Math.abs(value)) return -1;
           if(Math.abs(o.value) > Math.abs(value)) return 1;
           return 0;
       }


    }

}





