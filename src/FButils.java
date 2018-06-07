import java.util.ArrayList;
import java.util.Random;
import java.util.TreeMap;

public class FButils {
    public static final double eps = 1e-6;
    public static final double RANDOMTOLERANCE = 1e-6;
    public static double round(double d){
        return Math.floor( d + 0.5);
    }
    public static double randomRound(double d, double c){
        if(Math.abs(d - 0.5) <= RANDOMTOLERANCE){
            double p = Math.random();
            if(p > 0.5 && c>= 0) return 1d;
            if(p <= 0.5 && c>= 0) return 0d;
            if(p >= 0.5 ) return 1d;
            if(p < 0.5) return 0d;
        }
        return Math.floor( d + 0.5);
    }
    public static double randomRound(double d, double c, Random rand){
        if(Math.abs(d - 0.5) <= RANDOMTOLERANCE){
            double p = rand.nextDouble();
            if(p > 0.5 && c>= 0) return 1d;
            if(p <= 0.5 && c>= 0) return 0d;
            if(p >= 0.5 ) return 1d;
            if(p < 0.5) return 0d;
        }
        return Math.floor( d + 0.5);
    }

    public static double L2Norm(double [] d){
        double ret = 0;
        for(int i = 0 ; i < d.length ; i ++){
            ret = ret + d[i]*d[i];
        }
        ret = Math.sqrt(ret);
        return ret;
    }
    public static double L2Norm(ArrayList<Double> d){
        double ret = 0;
        for(Double di: d){
            ret = ret + di*di;
        }
        ret = Math.sqrt(ret);
        return ret;
    }

    public static double integralityGap(IntegralityGapTypes type, int niv, ArrayList<Double> x){
        double ret = 0d;

        if( type == IntegralityGapTypes.L2Norm){
            int i = 0;
            for( Double d : x){
                if(i >= niv){
                    break;
                }
                Double rd = FButils.round(d);
                ret = ret + Math.pow(Math.abs(rd - d),2);
                i++;
            }
            ret = Math.sqrt(ret)/(niv+0.0);
        }

        return ret;
    }
    public static double integralityGap(IntegralityGapTypes type, TreeMap<String,Double> x){
        double ret = 0d;

        if( type == IntegralityGapTypes.L2Norm){
            for( Double d : x.values()){
                Double rd = FButils.round(d);
                ret = ret + Math.pow(Math.abs(rd - d),2);
            }
            ret = Math.sqrt(ret)/(x.size()+0.0);
        }else if(type == IntegralityGapTypes.L1Norm){
            for( Double d : x.values()){
                Double rd = FButils.round(d);
                ret = ret + Math.abs(rd - d);
            }
            ret = ret/(x.size()+0.0);
        }

        return ret;
    }


    public static double objVal(ArrayList<Double> x, double[] c){
        double ret = 0;
        int i = 0;
        for(Double xi : x){
            ret = ret + xi*c[i];
            i++;
        }
        return ret;
    }

    public static double objVal(TreeMap<String, Double> xk, TreeMap<String, Double> xck,
                                TreeMap<String, Double> c, TreeMap<String, TreeMap<String, Double>> cq) {

        //IloLinearNumExpr expr = this.cplex.linearNumExpr();
        double objval = 0d;


        if(c != null && c.size() >= 1) {

            for(String xi : c.keySet()){
                if(xk.get(xi) != null) objval += c.get(xi)* xk.get(xi);
                else if(xck.get(xi) != null) objval += c.get(xi) * xck.get(xi);
            }
        }
        if(cq != null && cq.size() >= 1) {

            for(String xj : cq.keySet()) {
                for (String xi : cq.get(xj).keySet()) {
                    if (xck.get(xi) != null && xk.get(xj) != null)
                        objval += cq.get(xj).get(xi)* xck.get(xi)* xk.get(xj);
                    else if (xk.get(xi) != null && xk.get(xj) != null)
                        objval += cq.get(xj).get(xi)* xk.get(xi)* xk.get(xj);
                    else if (xck.get(xi) != null && xck.get(xj) != null)
                        objval += cq.get(xj).get(xi)*xck.get(xi)*xck.get(xj);
                    else if (xk.get(xi) != null && xck.get(xj) != null)
                        objval += cq.get(xj).get(xi)*xk.get(xi)*xck.get(xj);
                }
            }
        }

        return objval;
    }
    public static double objVal(TreeMap<String, Double> xk,
                                TreeMap<String, Double> c) {

        //IloLinearNumExpr expr = this.cplex.linearNumExpr();
        double objval = 0d;


        if(c != null && c.size() >= 1) {

            for(String xi : c.keySet()){
                if(xk.get(xi) != null) objval += c.get(xi)* xk.get(xi);
            }
        }


        return objval;
    }


    public static double L2norm(ArrayList<Double> xh, ArrayList<Double> x, int intvn){
        double ret = 0;
        for(int i = 0; i < intvn ; i++){
            ret = ret + Math.pow(x.get(i) - xh.get(i),2);
        }
        return ret;
    }
    public static double Chebynorm(ArrayList<Double> xh, ArrayList<Double> x, int intvn){
        double ret = 0;
        double val = 0;
        for(int i = 0; i < intvn ; i++){
            val = Math.abs(x.get(i) - xh.get(i));
            if(val > ret) ret = val;
        }
        return ret;
    }public static double L2norm(TreeMap<String, Double> xh, TreeMap<String,Double> x){
        double ret = 0d;
        for(String s : x.keySet()){
            ret = ret + Math.pow(x.get(s) - xh.get(s),2);
        }
        return ret;
    }
    public static double Chebynorm(TreeMap<String, Double> xh, TreeMap<String,Double> x){
        double ret = 0d;
        double val = 0d;
        for(String s : x.keySet()){
            val = Math.abs(x.get(s) - xh.get(s));
            if(val > ret) ret = val;
        }
        return ret;
    }
    public static double L1norm(ArrayList<Double> xh, ArrayList<Double> x, int intvn){
        double ret = 0;
        double val = 0;
        for(int i = 0; i < intvn ; i++){
            val = Math.abs(x.get(i) - xh.get(i));
             ret += val;
        }
        return ret;
    }
    public static double Chebynorm(ArrayList<Double> xh, int intvn){
        double ret = 0;
        double val = 0;
        for(int i = 0; i < intvn ; i++){
            val = Math.abs(xh.get(i));
            if(val > ret) ret = val;
        }
        return ret;
    }


    public static ArrayList<Double> sum(ArrayList<Double> xh, ArrayList<Double> x, int intvn){
        ArrayList<Double> ret = new ArrayList<Double>();
        double val = 0;
        for(int i = 0; i < intvn ; i++){
            val = (x.get(i) + xh.get(i));
            ret.add(val);
        }
        return ret;
    }

}
