import java.util.ArrayList;

public class FButils {
    public static final double eps = 1e-6;
    public static double round(double d){
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


    public static double objVal(ArrayList<Double> x, double[] c){
        double ret = 0;
        int i = 0;
        for(Double xi : x){
            ret = ret + xi*c[i];
            i++;
        }
        return ret;
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
