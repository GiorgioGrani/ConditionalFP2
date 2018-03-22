import ilog.concert.IloException;

import java.util.ArrayList;
//import java.util.List;

public class main {
    private static final double eps = 1e-8;
    public static void main(String [] aaa){
        int [] nsizes = {100};
        float [] pfactor = {5f};
        float [] perccon ={ 0.5f};
        int repeat = 1;

        for(int nsize : nsizes){
            for(float pf : pfactor){
                int p = Math.round(pf*nsize);
                for(float pc : perccon) {
                    for (int i = 0; i < repeat; i++) {
                        main.run(nsize, p);
                    }
                }
            }
        }
    }

    public static void run( int nsize , int p){
        ArrayList<Object> param = RANDOMINSIDEABOX( true, nsize, p);
        System.out.println("fine creazione");
        double [][] objectives =  (double[][]) param.get(0);
        double [][] matrixA = (double[][]) param.get(1);
        double [] b = (double[]) param.get(2);
        double [] lower = (double [] ) param.get(4);
        double [] upper = ( double [] ) param.get(5);

        boolean [] code = {true,false,false,false,false,true,false,false,true, true};
        if(code[0]){
            Simple fp = new Simple();
            try {
                fp.set(objectives[0], matrixA, b, lower, upper, nsize);
                ArrayList<Object> results = fp.solve();
                main.useResults(results, nsize, "Simple          ");
            }catch(IloException e){
                e.printStackTrace();
            }
        }
        if(code[1]) {
            LacosteJulien lj = new LacosteJulien();
            try {
                lj.set(objectives[0], matrixA, b, lower, upper, nsize);
                ArrayList<Object> results = lj.solve();
                main.useResults(results, nsize, "LacosteJulien          ");
            } catch (IloException e) {
                e.printStackTrace();
            }
        }
        if (code[2]) {
            LJ2 lj2 = new LJ2();
            try {
                lj2.set(objectives[0], matrixA, b, lower, upper, nsize);
                ArrayList<Object> results = lj2.solve();
                main.useResults(results, nsize, "OriginalLJ          ");
            } catch (IloException e) {
                e.printStackTrace();
            }
        }

        if(code[3]) {

            try {
                LJ3FS lj3 = new LJ3FS();
                lj3.set(objectives[0], matrixA, b, lower, upper, nsize);
                ArrayList<Object> results = lj3.solve();
                main.useResults(results, nsize, "NoArmijoLJ          ");
            } catch (IloException e) {
                e.printStackTrace();
            }
            System.gc();
        }
        if(code[4]) {
            try {
                LJ4S lj4 = new LJ4S();
                lj4.set(objectives[0], matrixA, b, lower, upper, nsize);
                ArrayList<Object> results = lj4.solve();
                main.useResults(results, nsize, "Sinusoidal          ");
            } catch (IloException e) {
                e.printStackTrace();
            }
        }
        if(code[5]) {
            try {
                BenchMark lj3 = new BenchMark();
                lj3.set(objectives[0], matrixA, b, lower, upper, nsize);
                ArrayList<Object> results = lj3.solve();
                main.useResults(results, nsize, "BenchMark          ");
            } catch (IloException e) {
                e.printStackTrace();
            }
        }
        if(code[6]) {
            try {
                FilledLJ5 lj3 = new FilledLJ5();
                lj3.set(objectives[0], matrixA, b, lower, upper, nsize);
                ArrayList<Object> results = lj3.solve();
                main.useResults(results, nsize, "Filled          ");
            } catch (IloException e) {
                e.printStackTrace();
            }

        }
        if(code[7]) {
            try {
                ConditionalLJ6 lj3 = new ConditionalLJ6();
                lj3.set(objectives[0], matrixA, b, lower, upper, nsize);
                ArrayList<Object> results = lj3.solve();
                main.useResults(results, nsize, "Conditional          ");
            } catch (IloException e) {
                e.printStackTrace();
            }

        }
        if(code[8]) {
            try {
                CLJ7 lj3 = new CLJ7();
                lj3.set(objectives[0], matrixA, b, lower, upper, nsize);
                ArrayList<Object> results = lj3.solve();
                main.useResults(results, nsize, "CLJ7          ");
            } catch (IloException e) {
                e.printStackTrace();
            }

        }
        if(code[9]) {
            try {
                CLJ8 lj3 = new CLJ8();
                lj3.set(objectives[0], matrixA, b, lower, upper, nsize);
                ArrayList<Object> results = lj3.solve();
                main.useResults(results, nsize, "CLJ8          ");
            } catch (IloException e) {
                e.printStackTrace();
            }

        }
    }

    public static ArrayList<Object> RANDOMINSIDEABOX(boolean allbinaries, int nsize, int pp){
        int m = 300;
        int n = nsize;//nsize;// (int) Math.round(Math.random()*200);
        int p = pp;//todo ricordati che hai fatto un knapsack 10*n;//(int) Math.round(Math.random()*200);
        //System.out.println(m+" "+n+" "+p);


        double [] [] matrix = new double [p] [n];
        double [] b = new double [p];
        double [][] o = new double [m][n];
        boolean [] binary = new boolean [n];
        double [] lower = new double[n];
        double [] upper = new double[n];

        for(int i = 0; i<p; i++){
            for(int j = 0; j<n; j++)
                matrix[i][j] = (Math.random()*21-1);
        }
        for(int i = 0; i<p; i++){
            b[i] =(10 - Math.random()*5);
            //System.out.println(b[i]);
        }

        double valup = 1;
        double vallow = 0;

        if(!allbinaries){
            valup = Math.round(1 * L1NormBounder(matrix, b));
            vallow = -valup;
        }

        for(int i =0 ; i< n; i++){
            lower[i] = vallow;
            upper[i] = valup;
        }
        //////////////////////////////////////





        for(int i = 0; i<m; i++){
            for(int j = 0; j<n; j++)
                o[i][j] = Math.round(Math.random()*100);//-50;
        }

        if(allbinaries){
            for(int i = 0 ; i < n ; i++){
                binary[i] = true;
            }
        }else {
            int count = 0;
            loop:
            while (true) {
                for (int i = 0; i < n; i++) {
                    if (Math.random() > 0.5d && !binary[i]) {
                        binary[i] = true;
                        count++;
                    }
                    if (count >= n / 2) {
                        break loop;
                    }

                }
            }
        }

        ArrayList<Object> ret = new ArrayList<>();
        ret.add(o);
        ret.add(matrix);
        ret.add(b);
        ret.add(binary);
        ret.add(lower);
        ret.add(upper);
        return ret;
    }


    public static double L1NormBounder( double[][] a, double [] b){
        double ret = 10e-200;
        for(int i = 0; i<a.length;i++){
            if(Math.abs(b[i]) > eps) {
                double maxb = b[i];
                for (int j = 0; j < a[0].length; j++) {
                    if (Math.abs(a[i][j]) > eps) {
                        double maxa = a[i][j];
                        double val = Math.abs(maxb/maxa);
                        if( val > ret){
                            ret = val;
                        }
                    }
                }
            }
        }
        return ret;
    }

    private static void useResults(ArrayList<Object> results, int nsize, String algorithm_name) {
        ArrayList<Double> x = (ArrayList<Double>) results.get(0);
        double ig = FButils.integralityGap(IntegralityGapTypes.L2Norm, nsize, x);
        boolean stop = (boolean) results.get(1);
        int iter = (int) results.get(2);
        long time = (long) results.get(3);
        boolean check = (boolean) results.get(4);
        double objgap = (double) results.get(5);
        System.out.println(algorithm_name+"  IntegralityGap: " + ig + "  Stop: " + stop + "  Iter: " + iter + "  Time (sec): " + (time / 1000) + "   Check:" + check + "  ObjGap: " + objgap);

    }
}
