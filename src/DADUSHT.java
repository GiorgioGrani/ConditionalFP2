import ilog.concert.*;
import ilog.cplex.IloCplex;
import java.util.Collections;
import java.util.Random;

import java.util.ArrayList;
import java.util.TreeMap;

public class DADUSHT  {
    private ArrayList<String> generals = new ArrayList<>();
    private ArrayList<String> binaries = new ArrayList<>();
    private ArrayList<String> continuous = new ArrayList<>();
    private TreeMap<String, IloNumVar> varbin = new TreeMap<String, IloNumVar>();
    private TreeMap<String, IloNumVar> vargen = new TreeMap<String, IloNumVar>();
    private TreeMap<String, IloNumVar> varcon = new TreeMap<String, IloNumVar>();
    private ArrayList<String> constraints = new ArrayList<>();
    private TreeMap<String, IloNumExpr> cons = new TreeMap<>();
    private TreeMap<String,  TreeMap <String, Double>> A = new TreeMap<>();
    private TreeMap<String,  TreeMap <String, Double>> V = new TreeMap<>();
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
    private double sign;
    private Random rand;

    private IloCplex cplex;
    private IloModel model;

    private IloAddable objective;
    private IloAddable distanceobj;
    private ArrayList<IloNumVar> x;


    public void set(ArrayList<Object > param ) throws IloException{
        this.rand = new Random(System.currentTimeMillis());


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
        this.sign = Math.pow(-1d, this.nbin + nint + 1);//*Math.pow(2,-this.numberOfIntegerVariables);

        this.createVariables();
        this.setObjective();
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
            //System.out.println(lb +" + " + ub);
        }
        for(String xi : generals){
            int lb = -(Integer.MAX_VALUE - 1);
            int ub =   (Integer.MAX_VALUE );
            if(lower.get(xi) != null) lb = (int) Math.round(lower.get(xi));
            if(upper.get(xi) != null) ub = (int)  Math.round(upper.get(xi));
            this.vargen.put(xi, this.cplex.numVar(lb,ub,xi));
        }
    }

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
        if(this.A == null || this.b == null) return;
        TreeMap<String, Double> norms = new TreeMap<>();
        for(String con : constraints){
            double norm = 0d;
            IloLQNumExpr constraint = this.cplex.lqNumExpr();


            if(this.A.get(con) != null && this.A.get(con).size() >= 1) {

                for(String xi : this.A.get(con).keySet()){
                    if(this.vargen.get(xi) != null)   constraint.addTerm(this.A.get(con).get(xi), this.vargen.get(xi));
                    else if(this.varbin.get(xi) != null) constraint.addTerm(this.A.get(con).get(xi), this.varbin.get(xi));
                    else if(this.varcon.get(xi) != null) constraint.addTerm(this.A.get(con).get(xi), this.varcon.get(xi));
                    norm += Math.pow(this.A.get(con).get(xi), 2);
                }
            }
            //todo ricordati che hai rimosso la parte quadratica


            if( directions.get(con) <= -0.5) this.cplex.addLe(constraint, this.b.get(con));
            else if( directions.get(con) >= 0.5) this.cplex.addGe(constraint, this.b.get(con));
            else this.cplex.addEq(constraint, this.b.get(con));
            this.cons.put(con, constraint);
            norms.put(con, norm);
        }

        double nmax = 1d;
        for(double d : norms.values()){
            if( nmax < d) nmax = d;
        }

        for(String con : this.constraints){
            for(String vn : this.A.get(con).keySet()){
            if(V.get(vn) == null) V.put(vn, new TreeMap<String, Double>());
                if(directions.get(con) > 0) {
                    V.get(vn).put(con,-A.get(con).get(vn)/(b.get(con)*nmax));
                }else if( directions.get(con) < 0){
                    V.get(vn).put(con,A.get(con).get(vn)/(b.get(con)*nmax));
                }else{
                    V.get(vn).put(con,A.get(con).get(vn)/(b.get(con)*nmax));
                }
            }
        }

    }


    private TreeMap<String, Double> GSwalk(ArrayList<String> A){
        Collections.sort(A);
        //double [] ret = new double [this.numberOfIntegerVariables];
        //double [][] B = new double [V.length][V[0].length];
        TreeMap<String, TreeMap<String, Double>> B = new TreeMap<>();
        ArrayList<String> subA = new ArrayList<>();
        for(String i : A){
            TreeMap<String, Double> pivot = new TreeMap<>();
            for(String s : V.get(i).keySet()){
                pivot.put(s,V.get(i).get(s));
            }
            B.put(i, pivot);

            for(String j : subA){
                double mu = 0;
                mu = scalarProd(B.get(j), B.get(i));
                double normb = Math.sqrt(scalarProd(B.get(j), B.get(j)));
                B.replace(i,update(B.get(i), mu/normb, B.get(j)));
            }
            subA.add(i);
        }

        String str = A.get(A.size() - 1);
        return B.get(str);
    }

    private TreeMap<String, Double> update(TreeMap<String, Double> vb, double mu, TreeMap<String, Double> vj){
        TreeMap<String, Double> ret = new TreeMap<>();
        for(String con : vb.keySet()){
            if(vj.get(con) != null) ret.put(con, vb.get(con) - mu * vj.get(con));
            else ret.put(con, vb.get(con));
        }
        for(String con : vj.keySet()){
            if(ret.get(con) == null) ret.put(con,  - mu * vj.get(con));
        }
        return ret;
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
            if(b.get(s) != null) ret += a.get(s)*b.get(s);
        }
        return ret;
    }

    private ArrayList<String> setA( TreeMap<String, Double> binaries){
        ArrayList<String> ret = new ArrayList<>();
        System.out.println("--------------- starterurururur");
        //System.out.println();
        for(String d : binaries.keySet()){
            //System.out.println("--> " + d);
            if(Math.abs(Math.round(binaries.get(d))-binaries.get(d)) > 0){
                System.out.println("--------------- "+d);
                ret.add(d);
            }
        }
        return ret;
    }

    private double [] deltaVal( TreeMap<String, Double> binaries,TreeMap<String, Double> weights, ArrayList<String> A){
        double dm = Double.MAX_VALUE -1;
        double dp =  (Double.MAX_VALUE -1 );


        for(String ind : A){
            double x = binaries.get(ind);
            //System.out.println(x+" <-"+ind);
            double val = weights.get(ind);
boolean check = true;
            if(Math.abs(val) > 0) {
                double d1 = -2d * x / val;
                double d2 = (2d - 2d * x) / val;

               for(String s : A){
                   double y = binaries.get(s);
                   System.out.println(ind+") d1 "+d1+" d2 " + d2+" x "+x +" y "+
                           y+" dd1 "+Math.abs(2d*y - 1d + d1*weights.get(s))+
                           " dd2 "+Math.abs(2d*y - 1d + d2*weights.get(s)));
                   if( Math.abs(2d*y - 1d + d1*weights.get(s)) >= (1 + eps) ||
                           Math.abs(2d*y - 1d + d2*weights.get(s)) >= (1 + eps)){
                       check = false;
                       break;
                   }
               }
               if (check){
                   double [] ret = new double [2];
                   ret[0] = d1*0.55;
                   ret[1] = d2*0.55;
                   return ret;
               }
            }
        }

        double [] ret = new double [2];
//        ret[0] = -dm;
//        ret[1] = dp;
        //Theoretical correct values
        ret[0] = -dm/2d;
        ret[1] = dp/2d;
System.out.println("ERROR FOUND");
        return ret;
    }

//    private double [] diff(double [] a , double [] b){
//        double [] ret = new double [a.length];
//        for(int i = 0; i < a.length; i++){
//            ret[i] = a[i] - b[i];
//        }
//        return ret;
//    }
private TreeMap<String, Double> diff(TreeMap<String, Double> vb,  TreeMap<String, Double> vj){
    TreeMap<String, Double> ret = new TreeMap<>();
    for(String con : vb.keySet()){
        if(vj.get(con) != null) ret.put(con, vb.get(con) -  vj.get(con));
        else ret.put(con, vb.get(con));
    }
    for(String con : vj.keySet()){
        if(ret.get(con) == null) ret.put(con,  -  vj.get(con));
    }
    return ret;
}

    private TreeMap<String, Double> solveSystem(TreeMap<String, Double> vnt,ArrayList<String> A) throws IloException{
        int n = A.size();
        //System.out.println(n);
        TreeMap<String, Double> rif = diff(vnt, V.get(A.get( n-1)));


        IloCplex sub = new IloCplex();
        sub.setOut(null);
        //sub.setParam(IloCplex.Param.Preprocessing.Presolve, false);
        TreeMap<String, IloNumVar> vars = new TreeMap<>();
        for(String s : A){
            if(! s.equalsIgnoreCase(A.get(n-1))) vars.put(s, sub.numVar(-(Double.MAX_VALUE-1), Double.MAX_VALUE));
        }

        for(String s : vnt.keySet()){
            IloLinearNumExpr expr = sub.linearNumExpr();
            for(String j : A){
                if(! j.equalsIgnoreCase(A.get(n-1))){
                    //System.out.println(n+" "+ind);
                    if(V.get(j).get(s) !=  null) expr.addTerm(V.get(j).get(s), vars.get(j));
                }
            }
            double rhs = 0d;
            if(rif.get(s) != null) rhs = rif.get(s);

            sub.addEq(expr, rhs);
        }


        IloLinearNumExpr obj = sub.linearNumExpr();
        for(IloNumVar v : vars.values()){
            obj.addTerm(0d, v); //todo puo' essere cambiato
        }
        sub.addMinimize(obj);
        sub.solve();

        TreeMap<String, Double>  ret = new TreeMap<>();

        for(String s : vars.keySet()){
            ret.put(s, sub.getValue(vars.get(s)));
        }
        ret.put(A.get(n-1), 1d);
        return ret;


    }

    private TreeMap<String, Double> update( TreeMap<String, Double> binaries, double [] delta,  TreeMap<String, Double> weights, ArrayList<String> A){
        double dm = delta[0];
        double dp = delta[1];
        double probm = dp/(dp + Math.abs(dm));
        //System.out.println(probm+"ffffffffffffffffffffffffffffffffffflllllllll");
        TreeMap<String, Double> xk = new TreeMap<>();

        double set = 0;
        for(String ind : A){
            if( rand.nextDouble() < probm){
                set = dm;
            }else{
                set = dp;
            }
            xk.put(ind, binaries.get(ind) + set+weights.get(ind) );
        }
        for(String s : binaries.keySet()){
            if(xk.get(s) == null) xk.put(s, binaries.get(s));
        }
        return xk;
    }

    private  TreeMap<String, Double> normalize(  TreeMap<String, Double> xk){
        double max = -(Double.MAX_VALUE - 1);
        double min = Double.MAX_VALUE;
        for(Double xi : xk.values()){
            if(xi > max) max = xi;
            if(xi < min) min = xi;
        }
        TreeMap<String, Double> newx = new  TreeMap<String, Double>();

        double c = (max - min) ;
        //System.out.println("C    "+c);
        //double rand = this.rand.nextDouble()*c + this.eps;
        //c = c + rand;
        //System.out.println(c);
        for(String  s : xk.keySet()){
            //newx.put(s,(xk.get(s) - min)/c);
            newx.put(s,(xk.get(s))/250d);
        }
        return newx;
    }

    public ArrayList<Object> solve() throws IloException{
        int maxiter = 0;
        this.cplex.setOut(null);
        this.cplex.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Barrier);
        //this.cplex.setParam(IloCplex.Param.Barrier.Limits.Iteration, 1);
        long start = System.currentTimeMillis();
        boolean solve = this.cplex.solve();
        //System.out.println(this.model);
        ArrayList<Object> ret = new ArrayList<>();
        ArrayList<Object> res = this.getX();
        TreeMap<String, Double> x =( TreeMap<String, Double>) res.get(0);
        TreeMap<String, Double> xc = ( TreeMap<String, Double>) res.get(1);
        TreeMap<String, Double> xk = x;
//        for(String s : x.keySet()){
//            x.replace(s,0.5);
//        }
     //   print(xk);

        boolean stop = this.stopCondition(xk);
        if(stop){
            long end = System.currentTimeMillis();
            ret.add(x);
            ret.add(stop);
            ret.add(0);
            ret.add(-start + end);
            ret.add(this.checkFeasibility(x, xc));
            ret.add( FButils.objVal(xk, this.c));
            return ret;
        }

        boolean print = false;

        TreeMap<String, Double> xh ;
        this.model.remove(this.objective);
//print(xk);
        ArrayList<String> A = setA(xk);

        while(!stop && maxiter < 100){
            //System.out.println(maxiter);
            maxiter++;
            TreeMap<String, Double> vnt = GSwalk(A);
            TreeMap<String, Double> weights;
            weights = solveSystem(vnt,A);
//            for(String s : weights.keySet()) System.out.print(s+" "+weights.get(s)+", ");
//            System.out.println();
            double [] delta ;
            delta = deltaVal(xk, weights, A);
            System.out.println("_____________________________start "+delta[0]+" "+delta[1]);

           print(xk);
            //System.out.println(delta[0]+" "+delta[1]+" "+maxiter);
            xk = update(xk, delta, weights, A);
//            System.out.println("___---___---0");
//            print(xk);
            //xk = normalize(xk);
//            System.out.println("___---___---1");
//            print(xk);
            xh = getXTilde(xk);
//            System.out.println("___---___---2");
//            print(xk);
//            System.out.println("___---___---2");
            if(checkFeasibility(xh, xc)){
                xk = xh;
                stop = true;
            }

            A = setA(xk);
            System.out.println(maxiter+"   ....   "+A.size());
            print(xk);

        }



        long end = System.currentTimeMillis();
        //print(xk);
        ret.add(xk);
        ret.add(stop);
        ret.add(maxiter);
        ret.add(-(start - end));
        ret.add(this.checkFeasibility(xk, xc));
        ret.add( FButils.objVal(xk, this.c));
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

