import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Collections;

public class BBFP  {
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

    private IloAddable setDistanceObjective( ArrayList<Double> Gx0) throws IloException{
        IloLinearNumExpr expr = this.cplex.linearNumExpr();

        int i = 0;
        for( Double d : Gx0){
            //double d = grad.get(i);
            expr.addTerm(d, this.x.get(i));
            i++;
            if(i >= this.numberOfIntegerVariables) break;
        }
        this.distanceobj = this.cplex.addMinimize(expr);
        return this.distanceobj;
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

    private IloAddable reinforcementCut(ArrayList<Double> G, ArrayList<Double> xh, double reinforce) throws IloException{
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
         return this.cplex.addGe(expr, reinforce);
        // );
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
        //System.out.println();
        for(int i = 0; i < this.numberOfIntegerVariables; i++){
            double d = g(xh.get(i))+g(xhm1.get(i));
            expr.addTerm(d, this.x.get(i));
        }
        //System.out.println();
        //System.out.println(this.cplex.addGe(expr, reinforce));

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
        if(indeces.size() == 0) return xh;
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


    private ArrayList<Object> branching(ArrayList<Double> xh, BranchingMode mode) throws IloException{
        int indmin = 0;

        if(mode == BranchingMode.maximumFranctionalValue) {
            double valmin = 0.5d;
            int i = 0;
            for (Double xi : xh) {
                double val = Math.abs(0.5d - xi);
                if (valmin > val) {
                    indmin = i;
                    valmin = val;
                } else if (val <= 1e-9) {
                    indmin = i;
                    valmin = val;
                    break;
                }
                if (i >= this.numberOfIntegerVariables-1) {
                    break;
                }
                i++;
            }
        }else if( mode == BranchingMode.A){
            int count = 0;
            int ind = 0;
            ArrayList<Integer> indeces = new ArrayList<>();
            for(Double xi : xh){
                if(Math.abs(xi-0.5) < 0.5-1e-10){
                    indeces.add(ind);
                    count++;
                }
                if (ind >= this.numberOfIntegerVariables-1) {
                    break;
                }
                ind++;
            }
            //ArrayList<Integer> consind = new ArrayList<>();
            int [] weights = new int[count];
            for(int i = 0 ; i < this.A.length; i ++){
                double val = 0;
                for( int j = 0 ; j < this.A[0].length; j++){
                    double xi = xh.get(j);
                    val += xi*A[i][j];
                }

                if(Math.abs(val - b[i]) <= 1e-11){
                    int base = 0;
                    for(Integer j : indeces){
                        if(Math.abs(A[i][j]) > 0){
                            weights[base] = weights[base] + 1;
                        }
                        base++;
                    }
                }
            }

            int maxval = 0;
            for(int i = 0 ; i < count; i++){
                if(maxval < weights[i]){
                    maxval = weights[i];
                    indmin = i;
                }
            }
            indmin = indeces.get(indmin);
        }else if( mode == BranchingMode.O){
            int count = 0;
            int ind = 0;
            ArrayList<Integer> indeces = new ArrayList<>();
            for(Double xi : xh){
                if(Math.abs(xi-0.5) < 0.5-1e-10){
                    indeces.add(ind);
                    count++;
                }

                if (ind >= this.numberOfIntegerVariables-1) {
                    break;
                }
                ind++;
            }
            //ArrayList<Integer> consind = new ArrayList<>();
            double [] weights = new double[count];
            for(int i = 0 ; i < this.A.length; i ++){
                double val = 0;
                int number = 0;

                for( int j = 0 ; j < this.A[0].length; j++){
                    double xi = xh.get(j);
                    val += xi*A[i][j];
                }
                for( Integer j : indeces){
                    if(Math.abs(A[i][j]) > 0){ number++;}
                }

                if(Math.abs(val - b[i]) <= 1e-11){
                    int base = 0;
                    for(Integer j : indeces){
                        if(Math.abs(A[i][j]) > 0){
                            weights[base] = weights[base] + A[i][j]/number;
                        }
                        base++;
                    }
                }
            }

            double maxval = 0;
            for(int i = 0 ; i < count; i++){
                if(maxval < weights[i]){
                    maxval = weights[i];
                    indmin = i;
                }
            }
            indmin = indeces.get(indmin);
        }
        //System.out.println("branching Variable :   "+indmin+"   vla:  "+this.cplex.getValue(this.x.get(indmin))+"  "+this.numberOfIntegerVariables);

        IloLinearNumExpr expr = this.cplex.linearNumExpr();
        expr.addTerm(1d, this.x.get(indmin));
        IloAddable [] res = new IloAddable [2];
        res[0] = this.cplex.addEq(expr, 0d);
        res[1] = this.cplex.addEq(expr, 1d);
        this.model.remove(res);
        ArrayList<Object> ret = new ArrayList<>();
        ret.add( res);
        ret.add( indmin);
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
        ArrayList<Double> poolpoint = new ArrayList<>();
        ArrayList<Double> olderpoint = new ArrayList<>();
        ArrayList<Double> G = G(xh);
        double valueobj = 1+scalarProdInt(G, xh);;
        this.model.remove(this.objective);
        int count = 0;

        LinkedList<ArrayList<IloAddable>> problems = new LinkedList<>();
        ArrayList<IloAddable> bunch = new ArrayList<>();
        LinkedList<IloAddable> objectives = new LinkedList<>();
        objectives.add(this.objective);
        IloAddable objective = this.objective;
        //int index = 0;
        //int actual = 0;
         int size = 1;
         int ites = 0;
        problems.addLast(bunch);
        boolean cast = true;

        while (!stop && maxiter <= 100 && problems.size() > 0) {
            maxiter++;
ites++;
            if(cast) {
                //bunch.add(
                 this.reinforcementCut(G, xh, valueobj);
                // );
                //this.model.remove(objective);
                //this.model.remove(this.distanceobj);
                objective = this.setDistanceObjective(G);
            }else{
                for(IloAddable con : bunch) {
                    this.model.remove(con);
                }
                this.model.remove(objective);
                //System.out.println(problems.size()+"   c");
                bunch = problems.getFirst();
                objective = objectives.getFirst();
                //System.out.println();
                for(IloAddable con : bunch) {
                    this.model.add(con);
                    //System.out.println(con);
                }
                this.model.add(objective);
                cast = true;
            }
            boolean solvable = this.cplex.solve();
            if(!solvable){
                problems.removeFirst();
                objectives.removeFirst();
                cast = false;
            }


            xk = this.getX();
            stop = this.stopCondition(xk);
            if(stop) break;

            //System.out.println("Here I am "+ count+" "+val+" "+(Math.abs(val) <= 1e-10)+" "+f(xk));
            xhm1 = xh;
            xh = this.getXTilde(xk);
            G = G(xh);
            boolean proximity = FButils.L2norm(xhm1, xh, numberOfIntegerVariables) <= eps;
            boolean inpool = false;
            if (!proximity){
                count = 0;
                valueobj = 1d + scalarProdInt(G, xh);
            } else if(proximity && count ==1) {
                cast = false;

                ArrayList<IloAddable> bunch1 = new ArrayList<>();
                ArrayList<IloAddable> bunch0 = new ArrayList<>();
                for(IloAddable con : bunch){
                    bunch1.add(con);
                    bunch0.add(con);
                }
                problems.removeFirst();
                objectives.removeFirst();
size += 2;
//System.out.println("FP iterations : "+(ites-1));
ites = -1;
                ArrayList<Object> results = branching(xk, BranchingMode.O);
                IloAddable [] cons = (IloAddable [] ) results.get(0);
                int index = (int) results.get(1);
                bunch0.add(cons[0]);
                bunch1.add(cons[1]);

                if(xh.get(index) >= 0.5){
                    problems.addFirst(bunch0);
                    problems.addLast(bunch1);
                    objectives.addFirst(objective);
                    objectives.addLast(objective);
                }else {
                    problems.addLast(bunch0);
                    problems.addFirst(bunch1);
                    objectives.addLast(objective);
                    objectives.addFirst(objective);
                }


            }else if (proximity  && count ==0) {
                for(IloAddable con : bunch){
                    this.model.remove(con);
                }
                this.cplex.solve();
                valueobj = Math.ceil(this.cplex.getObjValue());
                for(IloAddable con : bunch){
                    this.model.add(con);
                }
                count = 1;
            }

            boolean equals = FButils.L2norm(xhm1, xh, this.numberOfIntegerVariables) <= eps;


//todo se sono uguali non ha senso riproposse cplex
            if (this.n - this.numberOfIntegerVariables > 0) {
                this.model.remove(this.distanceobj);
                this.model.remove(objective);
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
                this.model.remove(objective);
            }


        }

        //System.out.println(size);
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

    private ArrayList<Double> getXTilde(ArrayList<Double> x) throws IloException{
        ArrayList<Double> ret = new ArrayList<>();
        boolean stop = true;
        for(int i = 0; i < this.numberOfIntegerVariables; i++){
            Double xi = x.get(i);
            Double rd = FButils.round(xi);
            //Double rd = FButils.randomRound(xi, this.c[i]);
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

