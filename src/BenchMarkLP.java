import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.util.List;
import java.util.TreeMap;
import java.util.Map;

import java.util.ArrayList;

public class BenchMarkLP  {

    private ArrayList<String> binaries = new ArrayList<>();
    private ArrayList<String> generals = new ArrayList<>();
    private ArrayList<String> continuous = new ArrayList<>();
    private TreeMap<String, IloNumVar> varbin = new TreeMap<String, IloNumVar>();
    private TreeMap<String, IloNumVar> vargen = new TreeMap<String, IloNumVar>();
    private TreeMap<String, IloNumVar> varcon = new TreeMap<String, IloNumVar>();
    private ArrayList<String> constraints = new ArrayList<>();
    private TreeMap<String,  TreeMap <String, Double>> A = new TreeMap<>();
    private TreeMap<String,  TreeMap <String,TreeMap<String,Double>>> Aq = new TreeMap<>();
    private TreeMap<String, Double> b = new TreeMap<>();
    private TreeMap<String, TreeMap<String, Double>> cq = new TreeMap<>();
    private TreeMap<String, Double> c = new TreeMap<>();
    private TreeMap<String, Double> lower = new TreeMap<>();
    private TreeMap<String, Double> upper = new TreeMap<>();
    private TreeMap<String, Integer> directions = new TreeMap<>();
    private TreeMap<String, IloNumExpr> cons = new TreeMap<>();

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
    //private IloAddable distanceobj;
    //private ArrayList< IloNumVar> x;

    private int tolerance;

    public void set(ArrayList<Object > param , int tolerance) throws IloException{
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
        this.setObjective();
        this.setConstraints();
        this.tolerance = tolerance;
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
            int lb =  0;
            int ub =  1;
            if(lower.get(xi) != null) lb =(int) Math.floor(0.5 + lower.get(xi));
            if(upper.get(xi) != null) ub =(int) Math.floor(0.5 + upper.get(xi));
            if( (ub-lb) < 1d-this.eps) this.varcon.put(xi, this.cplex.intVar(lb,ub,xi));
            else  this.varbin.put(xi, this.cplex.boolVar(xi));
        }
        for(String xi : generals){
            int lb = -(Integer.MAX_VALUE - 1);
            int ub =   (Integer.MAX_VALUE );
            if(lower.get(xi) != null) lb = (int) Math.round(lower.get(xi));
            if(upper.get(xi) != null) ub = (int)  Math.round(upper.get(xi));
            this.vargen.put(xi, this.cplex.intVar(lb,ub,xi));
        }
    }
    /*
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
*/


    //todo rimuovere objective dopo averlo utilizzato




    private void setObjective() throws IloException{

        //IloLinearNumExpr expr = this.cplex.linearNumExpr();
        IloLQNumExpr expr = this.cplex.lqNumExpr();


        if(this.c != null && this.c.size() >= 1) {

            for(String xi : this.c.keySet()){
                if(this.vargen.get(xi) != null)   expr.addTerm(this.c.get(xi), this.vargen.get(xi));
                else if(this.varbin.get(xi) != null) expr.addTerm(this.c.get(xi), this.varbin.get(xi));
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
//System.out.println(b.get(con)+" "+con);
//            for(String s : directions.keySet()) System.out.println(".."+s.hashCode());
            if( directions.get(con) <= -0.5) this.cplex.addLe(constraint, this.b.get(con));
            else if( directions.get(con) >= 0.5) this.cplex.addGe(constraint, this.b.get(con));
            else this.cplex.addEq(constraint, this.b.get(con));
            this.cons.put(con, constraint);
        }
    }









    private ArrayList<Object> getX() throws IloException{
        TreeMap<String, Double> x = new TreeMap<>();
        for(String b : varbin.keySet()){
            //System.out.println(b+" "+varbin.get(b));
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



   // private boolean checkFeasibility(ArrayList<Double> x){
//        if(this.A == null || this.b == null) return true;
//        for(int i = 0 ; i < this.A.length; i ++){
//            double val = 0;
//            for( int j = 0 ; j < this.A[0].length; j++){
//                double xi = x.get(j);
//                val += xi*A[i][j];
//            }
//
//            if(directions[i] > 0) {
//                if (val < b[i] - eps) return false;
//            }else if( directions[i] < 0){
//                if (val > b[i] + eps) return false;
//            }else{
//                if (val < b[i] - eps  || val > b[i] + eps  ) return false;
//            }
//        }
    //    return true;
    //}
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

           if( directions.get(con) <= -0.5 && (this.b.get(con) - conval) < -1e-9) return false;
           else if( directions.get(con) >= 0.5 && (this.b.get(con) - conval) > 1e-9) return false;
           else if( Math.abs(directions.get(con)) <= 0.5 &&  Math.abs(this.b.get(con) - conval) >= 1e-9)  return false;
       }
       return true;
   }





    public ArrayList<Object> solve() throws IloException{
        this.cplex.setOut(null);
        //this.cplex.setParam(IloCplex.Param.OptimalityTarget, 2);
        //this.cplex.setParam(IloCplex.DoubleParam.EpGap, 0d);
        //this.cplex.setParam(IloCplex.Param.MIP.Strategy.FPHeur, 1);
        //this.cplex.setParam(IloCplex.Param.Preprocessing.Presolve, false);
        if(this.tolerance >= 1) {
            this.cplex.setParam(IloCplex.Param.MIP.Limits.Solutions, this.tolerance);
        }

        long start = System.currentTimeMillis();
        //System.out.println(this.model);
        this.cplex.solve();
        //
        ArrayList<Object> ret = new ArrayList<>();
        TreeMap<String,Double> x = (TreeMap<String,Double> ) this.getX().get(0);
        TreeMap<String,Double> xc = (TreeMap<String,Double> ) this.getX().get(1);
        long end = System.currentTimeMillis();
        ret.add(x);
        ret.add(true);
        ret.add(0);
        ret.add(-start + end);
        ret.add(this.checkFeasibility(x,xc));
        ret.add(FButils.objVal(x,xc, this.c,this.cq));
        return ret;
    }

}
