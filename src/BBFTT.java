import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Collections;
import java.util.TreeMap;

public class BBFTT  {private ArrayList<String> binaries = new ArrayList<>();
    private ArrayList<String> generals = new ArrayList<>();
    private ArrayList<String> continuous = new ArrayList<>();
    private TreeMap<String, IloNumVar> varbin = new TreeMap<String, IloNumVar>();
    private TreeMap<String, IloNumVar> vargen = new TreeMap<String, IloNumVar>();
    private TreeMap<String, IloNumVar> varcon = new TreeMap<String, IloNumVar>();
    private ArrayList<String> constraints = new ArrayList<>();
    private TreeMap<String, IloNumExpr> cons = new TreeMap<>();
    private TreeMap<String,  TreeMap <String, Double>> A = new TreeMap<>();
    private TreeMap<String,  TreeMap <String,TreeMap<String,Double>>> Aq = new TreeMap<>();
    private TreeMap<String, Double> b = new TreeMap<>();
    private TreeMap<String, TreeMap<String, Double>> cq = new TreeMap<>();
    private TreeMap<String, Double> c = new TreeMap<>();
    private TreeMap<String, Double> lower = new TreeMap<>();
    private TreeMap<String, Double> upper = new TreeMap<>();
    private TreeMap<String, Integer> directions = new TreeMap<>();
    private String direction = "";
    private int n = 0;
    private int ncont = 0;
    private int nbin = 0;
    private int nint = 0;
    private int eq = 0;
    private int geq = 0;
    private int leq = 0;
    private static final double eps = 1e-6;

    private IloCplex cplex;
    private IloModel model;

    private IloAddable objective;
    private IloAddable distanceobj;

    public void set(ArrayList<Object > param ) throws IloException{
        this.n = ( int ) param.get(0);
        this.ncont = ( int ) param.get(1 );
        this.nbin = ( int ) param.get(2 );
        this.nint = ( int ) param.get(3 );
        this.eq = ( int ) param.get(4 );
        this.geq = ( int ) param.get(5 );
        this.leq = ( int ) param.get(6 );
        this.continuous = ( ArrayList<String> ) param.get(7 );
        this.binaries = ( ArrayList<String> ) param.get(8 );
        this.generals = ( ArrayList<String> ) param.get(9 );
        this.constraints = ( ArrayList<String> ) param.get(10 );
        this.direction = ( String ) param.get(11 );
        this.c = ( TreeMap<String, Double>  ) param.get(12 );
        this.cq = ( TreeMap<String, TreeMap<String, Double>> ) param.get(13 );
        this.A = ( TreeMap<String,  TreeMap <String, Double>> ) param.get(14 );
        this.Aq = ( TreeMap<String,  TreeMap <String,TreeMap<String,Double>>> ) param.get(15 );
        this.b = ( TreeMap<String, Double> ) param.get(16 );
        this.directions = ( TreeMap<String, Integer> ) param.get(17 );
        this.lower = ( TreeMap<String, Double> ) param.get(18 );
        this.upper = ( TreeMap<String, Double> ) param.get(19 );
        this.cplex = new IloCplex();
        this.model = this.cplex.getModel();

        this.createVariables();
        //this.setObjective();
        this.setConstraints();
        //System.out.println((model));
    }

    private void createVariables() throws IloException {
        for(String xi : continuous){
            double lb = - (Double.MAX_VALUE - 1d);
            double ub =   (Double.MAX_VALUE );
            if(lower.get(xi) != null) lb = lower.get(xi);
            if(upper.get(xi) != null) ub = upper.get(xi);
            this.varcon.put(xi, this.cplex.numVar(lb,ub,xi));
        }
        for(String xi : binaries){
            double lb =  0d;
            double ub =  1d;
            if(lower.get(xi) != null) lb = lower.get(xi);
            if(upper.get(xi) != null) ub = upper.get(xi);
            if( (ub-lb) < 1d-this.eps) this.varcon.put(xi, this.cplex.numVar(lb,ub,xi));
            else  this.varbin.put(xi, this.cplex.numVar(lb,ub,xi));
        }
        for(String xi : generals){
            int lb = -(Integer.MAX_VALUE - 1);
            int ub =   (Integer.MAX_VALUE );
            if(lower.get(xi) != null) lb = (int) Math.round(lower.get(xi));
            if(upper.get(xi) != null) ub = (int)  Math.round(upper.get(xi));
            this.vargen.put(xi, this.cplex.numVar(lb,ub,xi));
        }
    }


    //todo rimuovere objective dopo averlo utilizzato



    private void setObjective() throws IloException{

        //IloLinearNumExpr expr = this.cplex.linearNumExpr();
        IloLQNumExpr expr = this.cplex.lqNumExpr();
        //System.out.println(this.c);

        if(this.c != null && this.c.size() >= 1) {

            for(String xi : this.c.keySet()){
                if(this.vargen.get(xi) != null)   expr.addTerm(this.c.get(xi), this.vargen.get(xi));
                else if(this.varbin.get(xi) != null) {
                    //System.out.println("ciao mamma guarda come mi diverto");
                    expr.addTerm(this.c.get(xi), this.varbin.get(xi));
                }
                else if(this.varcon.get(xi) != null) expr.addTerm(this.c.get(xi), this.varcon.get(xi));
            }
        }
        if(this.cq != null && this.cq.size() >= 1) {

            for(String xj : this.cq.keySet()) {
                for (String xi : this.cq.get(xj).keySet()) {
                    if (this.vargen.get(xi) != null && this.vargen.get(xj) != null)
                        expr.addTerm(this.cq.get(xj).get(xi), this.vargen.get(xi), this.vargen.get(xj));
                    else if (this.vargen.get(xi) != null && this.varbin.get(xj) != null)
                        expr.addTerm(this.cq.get(xj).get(xi), this.vargen.get(xi), this.varbin.get(xj));
                    else if (this.vargen.get(xi) != null && this.varcon.get(xj) != null)
                        expr.addTerm(this.cq.get(xj).get(xi), this.vargen.get(xi), this.varcon.get(xj));
                    else if (this.varcon.get(xi) != null && this.varbin.get(xj) != null)
                        expr.addTerm(this.cq.get(xj).get(xi), this.varcon.get(xi), this.varbin.get(xj));
                    else if (this.varbin.get(xi) != null && this.varbin.get(xj) != null)
                        expr.addTerm(this.cq.get(xj).get(xi), this.varbin.get(xi), this.varbin.get(xj));
                    else if (this.varcon.get(xi) != null && this.varcon.get(xj) != null)
                        expr.addTerm(this.cq.get(xj).get(xi),this.varcon.get(xi),this.varcon.get(xj));
                    else if (this.varcon.get(xi) != null && this.vargen.get(xj) != null)
                        expr.addTerm(this.cq.get(xj).get(xi),this.varcon.get(xi),this.vargen.get(xj));
                    else if (this.varbin.get(xi) != null && this.varcon.get(xj) != null)
                        expr.addTerm(this.cq.get(xj).get(xi),this.varbin.get(xi),this.varcon.get(xj));
                    else if (this.varbin.get(xi) != null && this.vargen.get(xj) != null)
                        expr.addTerm(this.cq.get(xj).get(xi),this.varbin.get(xi),this.vargen.get(xj));
                }
            }
        }

        if(direction.equalsIgnoreCase("minimize"))
            this.objective = this.cplex.addMinimize(expr);
        else if(direction.equalsIgnoreCase("maximize"))
            this.objective = this.cplex.addMaximize(expr);
    }

    private void setConstraints() throws IloException{
        for(String con : constraints){
            IloLQNumExpr constraint = this.cplex.lqNumExpr();


            if(this.A.get(con) != null && this.A.get(con).size() >= 1) {

                for(String xi : this.A.get(con).keySet()){
                    if(this.vargen.get(xi) != null)   constraint.addTerm(this.A.get(con).get(xi), this.vargen.get(xi));
                    else if(this.varbin.get(xi) != null) constraint.addTerm(this.A.get(con).get(xi), this.varbin.get(xi));
                    else if(this.varcon.get(xi) != null) constraint.addTerm(this.A.get(con).get(xi), this.varcon.get(xi));
                }
            }
            if(this.Aq.get(con) != null && this.Aq.get(con).size() >= 1) {

                for(String xj : this.Aq.get(con).keySet()) {
                    for (String xi : this.Aq.get(con).get(xj).keySet()) {
                        if (this.vargen.get(xi) != null && this.vargen.get(xj) != null)
                            constraint.addTerm(this.Aq.get(con).get(xj).get(xi), this.vargen.get(xi), this.vargen.get(xj));
                        else if (this.vargen.get(xi) != null && this.varbin.get(xj) != null)
                            constraint.addTerm(this.Aq.get(con).get(xj).get(xi), this.vargen.get(xi), this.varbin.get(xj));
                        else if (this.vargen.get(xi) != null && this.varcon.get(xj) != null)
                            constraint.addTerm(this.Aq.get(con).get(xj).get(xi), this.vargen.get(xi), this.varcon.get(xj));
                        else if (this.varcon.get(xi) != null && this.varbin.get(xj) != null)
                            constraint.addTerm(this.Aq.get(con).get(xj).get(xi), this.varcon.get(xi), this.varbin.get(xj));
                        else if (this.varbin.get(xi) != null && this.varbin.get(xj) != null)
                            constraint.addTerm(this.Aq.get(con).get(xj).get(xi), this.varbin.get(xi), this.varbin.get(xj));
                        else if (this.varcon.get(xi) != null && this.varcon.get(xj) != null)
                            constraint.addTerm(this.Aq.get(con).get(xj).get(xi),this.varcon.get(xi),this.varcon.get(xj));
                        else if (this.varcon.get(xi) != null && this.vargen.get(xj) != null)
                            constraint.addTerm(this.Aq.get(con).get(xj).get(xi),this.varcon.get(xi),this.vargen.get(xj));
                        else if (this.varbin.get(xi) != null && this.varcon.get(xj) != null)
                            constraint.addTerm(this.Aq.get(con).get(xj).get(xi),this.varbin.get(xi),this.varcon.get(xj));
                        else if (this.varbin.get(xi) != null && this.vargen.get(xj) != null)
                            constraint.addTerm(this.Aq.get(con).get(xj).get(xi),this.varbin.get(xi),this.vargen.get(xj));
                    }
                }
            }

            if( directions.get(con) <= -0.5) this.cplex.addLe(constraint, this.b.get(con));
            else if( directions.get(con) >= 0.5) this.cplex.addGe(constraint, this.b.get(con));
            else this.cplex.addEq(constraint, this.b.get(con));
            this.cons.put(con, constraint);
        }
    }







    private ArrayList<Object> getX() throws IloException{
        TreeMap<String, Double> x = new TreeMap<>();
        for(String b : varbin.keySet()){
            x.put(b, this.cplex.getValue(this.varbin.get(b)));
        }
        TreeMap<String, Double> xc = new TreeMap<>();
        for(String b : varcon.keySet()){
            xc.put(b, this.cplex.getValue(this.varcon.get(b)));
        }
        ArrayList<Object> ret = new ArrayList<>();
        ret.add(x);
        ret.add(xc);
        return ret;
    }

    private double g(double d){
        if (Math.abs(d) <= this.eps) return 1d;
        if (Math.abs(1-d) <= this.eps) return -1d;
        return 0d;
    }

    private TreeMap<String,Double> G(TreeMap<String,Double> x0){
        TreeMap<String, Double> ret = new TreeMap<>();
        for(String d: x0.keySet()){
            ret.put(d,g(x0.get(d)));
        }
        return ret;
    }

    private IloAddable setDistanceObjective( TreeMap<String,Double> Gx0) throws IloException{
        IloLinearNumExpr expr = this.cplex.linearNumExpr();


        for( String d : Gx0.keySet()){
            //double d = grad.get(i);
            expr.addTerm(Gx0.get(d), this.varbin.get(d));
        }
        this.distanceobj = this.cplex.addMinimize(expr);
        return this.distanceobj;
    }

    private boolean checkFeasibility(TreeMap<String, Double> bin, TreeMap<String, Double> cont){

        for(String con : constraints){
            Double conval = 0d;


            if(this.A.get(con) != null && this.A.get(con).size() >= 1) {

                for(String xi : this.A.get(con).keySet()){
                    if(bin.get(xi) != null) conval += this.A.get(con).get(xi)*bin.get(xi);
                    else if(cont.get(xi) != null) conval += this.A.get(con).get(xi)*cont.get(xi);
                }
            }
            if(this.Aq.get(con) != null && this.Aq.get(con).size() >= 1) {

                for(String xj : this.Aq.get(con).keySet()) {
                    for (String xi : this.Aq.get(con).get(xj).keySet()) {
                        if (cont.get(xi) != null && bin.get(xj) != null)
                            conval += this.Aq.get(con).get(xj).get(xi)* cont.get(xi)* bin.get(xj);
                        else if (bin.get(xi) != null && bin.get(xj) != null)
                            conval += this.Aq.get(con).get(xj).get(xi)* bin.get(xi)* bin.get(xj);
                        else if (cont.get(xi) != null && cont.get(xj) != null)
                            conval += this.Aq.get(con).get(xj).get(xi)*cont.get(xi)*cont.get(xj);
                        else if (bin.get(xi) != null && cont.get(xj) != null)
                            conval += this.Aq.get(con).get(xj).get(xi)*bin.get(xi)*cont.get(xj);
                    }
                }
            }

            if( directions.get(con) <= -0.5 && (this.b.get(con) - conval) <= -1e-9) return false;
            else if( directions.get(con) >= 0.5 && (this.b.get(con) - conval) >= 1e-9) return false;
            else if( Math.abs(directions.get(con)) <= 0.5 &&  Math.abs(this.b.get(con) - conval) >= 1e-9)  return false;
        }
        return true;
    }



    private boolean stopCondition(TreeMap<String, Double> point){
        double val = 0d;
        for(Double d : point.values() ){ val += Math.pow( d - Math.floor(d + 0.5), 2);}
        if (val == Double.NaN) val = Double.MAX_VALUE;
        return val <= this.eps;
    }



    private double scalarProd( TreeMap<String, Double> a,  TreeMap<String, Double> b){

        double ret = 0d;
        for(String s : a.keySet()){
            ret += a.get(s)*b.get(s);
        }
        return ret;
    }

    private IloAddable reinforcementCut(TreeMap<String, Double> G,  double reinforce) throws IloException{
        IloLinearNumExpr expr = this.cplex.linearNumExpr();
        //System.out.println();
        for(String v : G.keySet()){
            expr.addTerm(G.get(v), this.varbin.get(v));
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


//
//
//
//    private void combinedReinforcementCut(ArrayList<Double> xh, double reinforce) throws IloException{
//        ArrayList<Object> check = checkProximality(xh);
//        boolean find = (boolean) check.get(0);
//        if(!find) return;
//        ArrayList<Double> xhm1 = (ArrayList<Double>) check.get(1);
//        IloLinearNumExpr expr = this.cplex.linearNumExpr();
//        //System.out.println();
//        for(int i = 0; i < this.numberOfIntegerVariables; i++){
//            double d = g(xh.get(i))+g(xhm1.get(i));
//            expr.addTerm(d, this.x.get(i));
//        }
//        //System.out.println();
//        //System.out.println(this.cplex.addGe(expr, reinforce));
//
//    }

//    private ArrayList<Double> localSearch(ArrayList<Double> xh){
//        ArrayList<Integer> indeces = new ArrayList<>();
//        int i = 0;
//        for(double xi : xh){
//            if( this.c[i] > 0 && Math.abs(xi) <= this.eps){
//                indeces.add(i);
//            }else if( this.c[i] < 0 && Math.abs(xi) >= 1d-this.eps){
//                indeces.add(i);
//            }
//            i++;
//            if (i >= this.numberOfIntegerVariables) break;
//        }
//        if(indeces.size() == 0) return xh;
//        int minimumind = indeces.get(0);
//        for(int ind : indeces){
//            if(Math.abs(this.c[minimumind]) < Math.abs(this.c[ind])) minimumind = ind;
//        }
//
//        ArrayList<Double> ret = new ArrayList<>();
//        int j = 0;
//        for(Double xi : xh){
//            double add = xi;
//            if(j == minimumind){
//                if(Math.abs(xi) >= 1d-this.eps) ret.add(0d);
//                else ret.add(1d);
//            }else{
//                ret.add(xi);
//            }
//            j++;
//            if(j >= this.numberOfIntegerVariables) break;
//        }
//        return ret;
//
//    }


    private ArrayList<Object> branching(TreeMap<String, Double> xh, BranchingMode mode) throws IloException{
        //int indmin = 0;
        String sret = "";

        if(mode == BranchingMode.maximumFranctionalValue) {
            double valmin = 0.5d;
            //int i = 0;
            for (String s : xh.keySet()) {
                double val = Math.abs(0.5d - xh.get(s));
                if (valmin > val) {
                    sret = s;
                    valmin = val;
                } else if (val <= 1e-9) {
                    sret = s;
                    valmin = val;
                    break;
                }
            }
        }else if( mode == BranchingMode.A){
            int count = 0;
            //int ind = 0;
            ArrayList<String> indeces = new ArrayList<>();
            for(String s : xh.keySet()){
                if(Math.abs(xh.get(s)-0.5) < 0.5-1e-10){
                    indeces.add(s);
                    count++;
                }
            }
            //ArrayList<Integer> consind = new ArrayList<>();

            TreeMap<String , Integer> weights = new TreeMap<>();
            for(String s : indeces) weights.put(s, 0);
            for(String name : this.constraints){
                double val = this.cplex.getValue(this.cons.get(name));
                double b = this.b.get(name);
                if( Math.abs(val - b) <= 1e-9){
                    for(String s : indeces){
                        if( (this.A.get(name) != null && this.A.get(name).get(s) != null && Math.abs(this.A.get(name).get(s)) >= 1e-9)) weights.replace(s, weights.get(s) + 1);
                        else if( (this.Aq.get(name) != null && this.Aq.get(name).get(s) != null) ){
                            for(String v : this.Aq.get(name).get(s).keySet()) {
                                if( Math.abs(  this.Aq.get(name).get(s).get(v)) >= 1e-9 ) {
                                    weights.replace(s, weights.get(s) + 1);
                                }
                            }
                        }else if (this.Aq.get(name) != null){
                            for(String z : this.Aq.get(name).keySet()){
                                for(String v : this.Aq.get(name).get(z).keySet()) {
                                    if( z.equalsIgnoreCase(v) &&
                                            Math.abs(  this.Aq.get(name).get(z).get(v)) >= 1e-9 ){
                                        weights.replace(s, weights.get(s) + 1);
                                    }
                                }
                            }
                        }


                    }
                }
            }

            int val = 0;
            String guard = "";
            for(String s : weights.keySet()){
                guard = s;
                if(weights.get(s) > val){
                    val = weights.get(s);
                    sret = s;
                }
            }
            if(sret.equals("")){
                sret = guard;
            }
        }
//        else if( mode == BranchingMode.O){
//            int count = 0;
//            int ind = 0;
//            ArrayList<Integer> indeces = new ArrayList<>();
//            for(Double xi : xh){
//                if(Math.abs(xi-0.5) < 0.5-1e-10){
//                    indeces.add(ind);
//                    count++;
//                }
//
//
//            }
//            //ArrayList<Integer> consind = new ArrayList<>();
//            double [] weights = new double[count];
//            for(int i = 0 ; i < this.A.length; i ++){
//                double val = 0;
//                int number = 0;
//
//                for( int j = 0 ; j < this.A[0].length; j++){
//                    double xi = xh.get(j);
//                    val += xi*A[i][j];
//                }
//                for( Integer j : indeces){
//                    if(Math.abs(A[i][j]) > 0){ number++;}
//                }
//
//                if(Math.abs(val - b[i]) <= 1e-11){
//                    int base = 0;
//                    for(Integer j : indeces){
//                        if(Math.abs(A[i][j]) > 0){
//                            weights[base] = weights[base] + A[i][j]/number;
//                        }
//                        base++;
//                    }
//                }
//            }
//
//            double maxval = 0;
//            for(int i = 0 ; i < count; i++){
//                if(maxval < weights[i]){
//                    maxval = weights[i];
//                    indmin = i;
//                }
//            }
//            indmin = indeces.get(indmin);
//        } todo reimplementare O
        //System.out.println("branching Variable :   "+indmin+"   vla:  "+this.cplex.getValue(this.x.get(indmin))+"  "+this.numberOfIntegerVariables);

        IloLinearNumExpr expr = this.cplex.linearNumExpr();
        expr.addTerm(1d, this.varbin.get(sret));
        IloAddable [] res = new IloAddable [2];
        res[0] = this.cplex.addEq(expr, 0d);
        res[1] = this.cplex.addEq(expr, 1d);
        this.model.remove(res);
        ArrayList<Object> ret = new ArrayList<>();
        ret.add( res);
        ret.add( sret);
        return ret;

    }

    public ArrayList<Object> solve() throws IloException {
        int maxiter = 0;
        this.cplex.setOut(null);
        this.cplex.setParam(IloCplex.Param.OptimalityTarget, 2);
        //this.cplex.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Barrier);

        long start = System.currentTimeMillis();
        //System.out.println(model);
        //this.setObjective();
        this.cplex.solve();

        //System.out.println("giggino   " + cplex.getObjective());
        //System.out.println(this.model);
        ArrayList<Object> ret = new ArrayList<>();

        ArrayList<Object> get = this.getX();
        TreeMap<String, Double> x = (TreeMap<String, Double>) get.get(0);
        TreeMap<String, Double> xc = (TreeMap<String, Double>) get.get(1);
        TreeMap<String, Double> xk = x;
        TreeMap<String, Double> xck = xc;


        boolean stop = this.stopCondition(xk);
        boolean feas = false;
        double lpobj = this.cplex.getObjValue();
        if (stop) {
            //System.out.println("non mi sono fermato"); todo trovare perché torna subito stop true anche se l'integralità non è soddisfatta
            long end = System.currentTimeMillis();
            ret.add(x);
            ret.add(stop);
            ret.add(-1);
            ret.add(-start + end);
            ret.add(this.checkFeasibility(x,xc));
            ret.add(FButils.objVal(xk,xck, this.c, this.cq));
            //print(x);
            return ret;
        }

        this.setObjective();
        this.cplex.solve();

        get = this.getX();
        x = (TreeMap<String, Double>) get.get(0);
        xc = (TreeMap<String, Double>) get.get(1);
        xk = x;
        xck = xc;


        stop = this.stopCondition(xk);
        feas = false;
        lpobj = this.cplex.getObjValue();
        if (stop) {
            //System.out.println("non mi sono fermato"); todo trovare perché torna subito stop true anche se l'integralità non è soddisfatta
            long end = System.currentTimeMillis();
            ret.add(x);
            ret.add(stop);
            ret.add(-1);
            ret.add(-start + end);
            ret.add(this.checkFeasibility(x,xc));
            ret.add(FButils.objVal(xk,xck, this.c, this.cq));
            //print(x);
            return ret;
        }

        boolean print = false;

        TreeMap<String, Double> xh = this.getXTilde(xk);
        TreeMap<String, Double> xhm1 = new TreeMap<>();
        //ArrayList<Double> poolpoint = new ArrayList<>();
        //ArrayList<Double> olderpoint = new ArrayList<>();
        TreeMap<String, Double> G = G(xh);
        double valueobj = 1 + scalarProd(G, xh);;
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

        while (!stop && maxiter <= 1000 && problems.size() > 0) {
            maxiter++;
            ites++;
            if(cast) {
                //bunch.add(
                this.reinforcementCut(G, valueobj);
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
                continue;
            }


            get = this.getX();
            xk = (TreeMap<String, Double>) get.get(0);
            xck = (TreeMap<String, Double>) get.get(1);
            stop = this.stopCondition(xk);
            if(stop) break;

            //System.out.println("Here I am "+ count+" "+val+" "+(Math.abs(val) <= 1e-10)+" "+f(xk));
            xhm1 = xh;
            xh = this.getXTilde(xk);
            G = G(xh);
            boolean proximity = FButils.L2norm(xhm1, xh) <= eps;
            //boolean inpool = false;
            if (!proximity){
                count = 0;
                valueobj = 1d + scalarProd(G, xh);
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
                ArrayList<Object> results = branching(xk, BranchingMode.A);
                IloAddable [] cons = (IloAddable [] ) results.get(0);
                String index = (String) results.get(1);
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

            //boolean equals = FButils.L2norm(xhm1, xh, this.numberOfIntegerVariables) <= eps;


//todo se sono uguali non ha senso riproposse cplex
            if (this.ncont > 0) {
                this.model.remove(this.distanceobj);
                this.model.remove(objective);
                this.model.add(this.objective);
                for (String s : xh.keySet()) {
                    double value = xh.get(s);
                    this.varbin.get(s).setLB(value);
                    this.varbin.get(s).setUB(value);
                }
                feas = this.cplex.solve();
                //System.out.println("feas  " + feas);
                if (feas) {
                    for (String s : xck.keySet()) {
                        xck.replace(s,this.cplex.getValue(this.varcon.get(s)));
                    }
                    xk = xh;
                    stop = true;
                    break;
                }


                for (String s : xh.keySet()) {
                    this.varbin.get(s).setLB((this.lower.get(s) != null ? this.lower.get(s) : 0d));
                    this.varbin.get(s).setUB((this.upper.get(s) != null ? this.upper.get(s) : 1d));
                }
                this.model.remove(this.objective);
            } else if (checkFeasibility(xh,xck)) {
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
        ret.add(this.checkFeasibility(xk,xck));
        ret.add(FButils.objVal(xk,xck, this.c,this.cq));
        return ret;
    }

    private TreeMap<String, Double> getXTilde(TreeMap<String, Double> x) throws IloException{
        TreeMap<String, Double> ret = new TreeMap<>();
        boolean stop = true;
        for(String s : x.keySet()){
            Double xi = x.get(s);
            Double rd = FButils.round(xi);
            //Double rd = FButils.randomRound(xi, this.c[i]);
            ret.put(s,rd);
            //System.out.println(rd+"   "+xi);
        }
        return ret;
    }

    private void print(TreeMap<String, Double> xk){
        for(String s : xk.keySet()){
            System.out.println(s+"                         "+xk.get(s));
        }
    }

}

