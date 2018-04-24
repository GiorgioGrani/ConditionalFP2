import com.csvreader.CsvWriter;
import ilog.concert.IloException;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import  java.util.Scanner;
import  java.io.File;
import java.io.FileNotFoundException;
import java.util.TreeMap;
//import java.util.List;

public class main {
    private static final double eps = 1e-8;
    public static void main(String [] aaa){

//        int [] nsizes = {1000};
//        float [] pfactor = {2f};
//        float [] perccon ={ 0.0f};
//        int repeat = 1;
//
//        for(int nsize : nsizes){
//            for(float pf : pfactor){
//                int p = Math.round(pf*nsize);
//                for(float pc : perccon) {
//                    int intvar = Math.round((1f-pc)*nsize);
//                    int convar = Math.max( 0, nsize - intvar);
//                    for (int i = 0; i < repeat; i++) {
//                        main.run(intvar, convar, p);
//                    }
//                }
//            }
//        }

        run( aaa[0], aaa[1]);
    }

    public static ArrayList<Object> readMPS(String mps) throws FileNotFoundException{
        Scanner scan = new Scanner(new File(mps));
        String name = mps;
        TreeMap<String, String> con_name_type = new TreeMap<>();
        String next = scan.next();
        if(next.equalsIgnoreCase("NAME")){
            name = scan.next();
            next = scan.next();
        }
        if(next.equalsIgnoreCase("ROWS")){
            next = scan.next();
            while( !next.equalsIgnoreCase("COLUMNS")){
                String direction = next;
                String conname = scan.next();
                con_name_type.put(conname, direction);
                next = scan.next();
            }
        }
        TreeMap<String, TreeMap<String, Double>> intvars = new TreeMap<>();
        TreeMap<String, TreeMap<String, Double>> numvars = new TreeMap<>();

        if (next.equalsIgnoreCase("COLUMNS")) {
            next = scan.next();
            while (!next.equalsIgnoreCase("RHS") && !next.equalsIgnoreCase("BOUNDS")) {
                if (next.equalsIgnoreCase("MARK0000") || next.equalsIgnoreCase("MARK") ) {
                    scan.next();
                    next = scan.next();
                    next = scan.next();
                    while (!next.equalsIgnoreCase("MARK0001") && !next.equalsIgnoreCase("MARKEND")) {


                        String var = next;
                        TreeMap<String, Double> submap = new TreeMap<>();

                        while (true) {
                            next = scan.next();
                            //System.out.println(next);
                            if (scan.hasNextDouble()) {
                                Double d = scan.nextDouble();
                                submap.put(next, d);
                            } else if (next.equalsIgnoreCase(var)) {

                            } else {
                                break;
                            }
                        }
                        intvars.put(var, submap);

                    }

                    scan.next();
                    scan.next();
                    next = scan.next();
                }
                if (!next.equalsIgnoreCase("RHS")) {
                    while (!next.equalsIgnoreCase("RHS") && !next.equalsIgnoreCase("BOUNDS") && !next.equalsIgnoreCase("MARK0000") && !next.equalsIgnoreCase("MARK")) {


                        String var = next;
                        TreeMap<String, Double> submap = new TreeMap<>();

                        while (true) {
                            next = scan.next();
                            if (scan.hasNextDouble()) {
                                Double d = scan.nextDouble();
                                submap.put(next, d);
                                //System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@           "+d);
                            } else if (next.equalsIgnoreCase(var)) {

                            } else {
                                break;
                            }
                        }
                        numvars.put(var, submap);

                    }

                }

            }

        }

//        if (!next.equalsIgnoreCase("RHS")) {
//            while (!next.equalsIgnoreCase("RHS") && !next.equalsIgnoreCase("BOUNDS")) {
//
//
//                String var = next;
//                TreeMap<String, Double> submap = new TreeMap<>();
//
//                while (true) {
//                    next = scan.next();
//                    if (scan.hasNextDouble()) {
//                        Double d = scan.nextDouble();
//                        submap.put(next, d);
//                        //System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@           "+d);
//                    } else if (next.equalsIgnoreCase(var)) {
//
//                    } else {
//                        break;
//                    }
//                }
//                numvars.put(var, submap);
//
//            }
//
//        }


        TreeMap<String, Double> rhs = new TreeMap<>();
        if (next.equalsIgnoreCase("RHS")) {
            next = scan.next();
            while (next.equalsIgnoreCase("rhs")) {
                //System.out.println("                                        )) " +next);
                rhs.put(scan.next(), scan.nextDouble());
                next = scan.next();
                if(scan.hasNextDouble()){
                    rhs.put(next, scan.nextDouble());
                    next = scan.next();
                }
                //
            }
        }
        TreeMap<String, Double> upperbound = new TreeMap<>();
        TreeMap<String, Double> lowerbound = new TreeMap<>();

        if (next.equalsIgnoreCase("BOUNDS")) {
            //next = scan.next();
            while (scan.hasNext()) {
                next = scan.next();
                //System.out.println(next);
                if (next.equalsIgnoreCase("UP") || next.equalsIgnoreCase("UI")) {
                    scan.next();
                    String xname = scan.next();
                    upperbound.put(xname, scan.nextDouble());
                   // System.out.println("upupupupupup   "+xname);
                } else if (next.equalsIgnoreCase("LO") || next.equalsIgnoreCase("LI")) {

                    scan.next();
                    String xname = scan.next();
                    lowerbound.put(xname, scan.nextDouble());
                    //System.out.println("ooooomionodno   "+xname);
                } else if (next.equalsIgnoreCase("FX")) {
                    scan.next();
                    String xname = scan.next();
                    Double xval = scan.nextDouble();
                    lowerbound.put(xname, xval);
                    upperbound.put(xname, xval);
                    //System.out.println(xname+" yyyyyyyyyyyyyyyyyyyyyy  "+xval);
                } else if (next.equals("FR")) {
                    scan.next();
                    String xname = scan.next();
                    //Double xval = scan.nextDouble();
                    lowerbound.put(xname, -(Double.MAX_VALUE-1));
                    upperbound.put(xname, Double.MAX_VALUE);
                    //System.out.println(xname+" yyyyyyyyyyyyyyyyyyyyyy  "+xval);
                } else if (next.equalsIgnoreCase("BV")) {
                    scan.next();
                    String xname = scan.next();
                    //Double xval = scan.nextDouble();
                    lowerbound.put(xname, 0d);
                    upperbound.put(xname, 1d);
                    //System.out.println(xname+" yyyyyyyyyyyyyyyyyyyyyy  "+xval);
                } else if (next.equalsIgnoreCase("MI")) {
                    scan.next();
                    String xname = scan.next();
                    //Double xval = scan.nextDouble();
                    lowerbound.put(xname, -(Double.MAX_VALUE-1));
                    upperbound.put(xname, 0d);
                    //System.out.println(xname+" yyyyyyyyyyyyyyyyyyyyyy  "+xval);
                }else if (next.equalsIgnoreCase("PL")) {
                    scan.next();
                    String xname = scan.next();
                    //Double xval = scan.nextDouble();
                    lowerbound.put(xname, 0d);
                    upperbound.put(xname, Double.MAX_VALUE);
                    //System.out.println(xname+" yyyyyyyyyyyyyyyyyyyyyy  "+xval);
                }else if(next.equalsIgnoreCase("ENDATA")){
                    break;
                }
            }
        }


        TreeMap<String, Integer> codvars = new TreeMap<>();
        TreeMap<String, Integer> codcons = new TreeMap<>();

        int i = 0;
        for(String xi : intvars.keySet()){
            codvars.put(xi, i);
            i++;
        }
        for(String xi : numvars.keySet()){
            codvars.put(xi, i);
            i++;
        }


        int objval = 0;
        i = 0;
        for(String xi : con_name_type.keySet()){
            codcons.put(xi, i);
            if(con_name_type.get(xi).equalsIgnoreCase("N")) objval = i;
            i++;
        }


        int numintvar = intvars.size();
        int n = numintvar + numvars.size();
        int p = con_name_type.size() - 1;

        double [] [] matrix = new double [p] [n];
        double [] b = new double [p];
        int [] directions = new int [p];
        double [] o = new double [n];
        boolean [] binary = new boolean [n];
        double [] lower = new double[n];
        double [] upper = new double[n];
        for(int r = 0; r< upper.length; r++){
            upper[r] = Double.MAX_VALUE;
        }

        for(String xi : intvars.keySet()){
            int k = codvars.get(xi);
            for(String ci : intvars.get(xi).keySet()){
                int h = 0;
                double val = 0;
                for(String con : codcons.keySet()){

                    if(con.equalsIgnoreCase(ci)){
                        h = codcons.get(con);
                        break;
                    }
                }
                val = intvars.get(xi).get(ci);
                if(h==objval){
                    o[k] = val;
                }else {
                    if( h > objval) matrix[h-1][k] = val;
                    else matrix[h][k] = val;

                }
            }
        }
        for(String xi : numvars.keySet()){
            int k = codvars.get(xi);
            for(String ci : numvars.get(xi).keySet()){
                int h = 0;
                double val = 0;
                for(String con : codcons.keySet()){

                    if(con.equalsIgnoreCase(ci)){
                        h = codcons.get(con);
                        break;
                    }
                }
                val = numvars.get(xi).get(ci);
                if(h==objval){
                    o[k] = val;
                }else {
                    if( h > objval) matrix[h-1][k] = val;
                    else matrix[h][k] = val;

                }
            }
        }

        for(String ci : rhs.keySet()){
            int h = 0;
            for(String s : codcons.keySet()){
                if(s.equalsIgnoreCase(ci)){
                    h = codcons.get(s);
                    break;
                }
            }
            if (h != objval) {
                if (h > objval) b[h - 1] = rhs.get(ci);
                else b[h] = rhs.get(ci);
            }

        }
        for (String ci : con_name_type.keySet()) {
            int h = 0;
            for (String s : codcons.keySet()) {
                if (s.equalsIgnoreCase(ci)) {
                    h = codcons.get(s);
                    break;
                }
            }
            if (h != objval) {
                if (h > objval) directions[h - 1] = codify(con_name_type.get(ci));
                else directions[h] = codify(con_name_type.get(ci));
            }
        }

        for(String xi : upperbound.keySet()){
            int k = 0;
            for(String s : codvars.keySet()){
                if(s.equalsIgnoreCase(xi)){
                    k = codvars.get(s);
                    break;
                }
            }

            upper[k] = upperbound.get(xi);

        }

        for(String xi : lowerbound.keySet()){
            int k = 0;
            for(String s : codvars.keySet()){
                if(s.equalsIgnoreCase(xi)){
                    k = codvars.get(s);
                    break;
                }
            }

            lower[k] = lowerbound.get(xi);

        }
        

        ArrayList<Object> param = new ArrayList<>();
        param.add(numintvar);
        param.add(o);
        param.add(matrix);
        param.add(b);
        param.add(lower);
        param.add(upper);
        param.add(directions);
        return param;
    }

    public static int codify(String s){
        if(s.equalsIgnoreCase("L")) return -1;
        else if(s.equalsIgnoreCase("E")) return 0;
        else return 1;
    }

    public static void run( String mps, String output){

        ArrayList<Object> param = new ArrayList<>();
        try {
            param = readMPS(mps);
        }catch (FileNotFoundException e){
            e.printStackTrace();
            return;
        }
        int nsize = (int) param.get(0);
        System.out.println("fine creazione");
        double [] objectives =  (double []) param.get(1);
        double [][] matrixA = (double[][]) param.get(2);
        double [] b = (double[]) param.get(3);
        double [] lower = (double [] ) param.get(4);
        double [] upper = ( double [] ) param.get(5);
        int [] directions = ( int [] ) param.get(6);

//        printMatrix(objectives);
//        printMatrix(matrixA);
//        printMatrix(b);
//        printMatrix(lower);
//        printMatrix(upper);
//        printMatrix(directions);

        ArrayList<Double> xb = new ArrayList<>();
        ArrayList<Double> xc = new ArrayList<>();

        boolean [] code = {true,true,true,true,false,true,true,true,false , false, false};
        if(code[0]){
            BasicVersion fp = new BasicVersion();
            try {
                fp.set(objectives, matrixA, b, lower, upper, nsize, directions);
                ArrayList<Object> results = fp.solve();
                main.useResults(results, nsize, "Simple", output);
            }catch(IloException e){
                e.printStackTrace();
            }
        }
        if(code[1]) {
            ReinforcementFP fp = new ReinforcementFP();
            try {
                fp.set(objectives, matrixA, b, lower, upper, nsize, directions);
                ArrayList<Object> results = fp.solve();
                main.useResults(results, nsize, "Reinforcement", output);
            }catch(IloException e){
                e.printStackTrace();
            }
        }
        if(code[2]) {
            AggressiveRFP fp = new AggressiveRFP();
            try {
                fp.set(objectives, matrixA, b, lower, upper, nsize, directions);
                ArrayList<Object> results = fp.solve();
                main.useResults(results, nsize, "Aggressive", output);
            }catch(IloException e){
                e.printStackTrace();
            }
        }

        if(code[3]) {

            CDFP fp = new CDFP();
            try {
                fp.set(objectives, matrixA, b, lower, upper, nsize, directions);
                ArrayList<Object> results = fp.solve();
                main.useResults(results, nsize, "CDFP", output);
            }catch(IloException e){
                e.printStackTrace();
            }
        }
        if(code[4]) {
            try {
                LJ4S lj4 = new LJ4S();
                lj4.set(objectives, matrixA, b, lower, upper, nsize);
                ArrayList<Object> results = lj4.solve();
                //main.useResults(results, nsize, "Sinusoidal          ");
            } catch (IloException e) {
                e.printStackTrace();
            }
        }
        if(code[5]) {
            try {
                BenchMark lj3 = new BenchMark();
                lj3.set(objectives, matrixA, b, lower, upper, nsize, directions, 2);
                ArrayList<Object> results = lj3.solve();
                main.useResults(results, nsize, "BenchMark2", output);
                xb = (ArrayList<Double>) results.get(0);
            } catch (IloException e) {
                e.printStackTrace();
            }
        }
        if(code[6]) {
            try {
                BenchMark lj3 = new BenchMark();
                lj3.set(objectives, matrixA, b, lower, upper, nsize, directions,3 );
                ArrayList<Object> results = lj3.solve();
                main.useResults(results, nsize, "BenchMark3", output);
                xb = (ArrayList<Double>) results.get(0);
            } catch (IloException e) {
                e.printStackTrace();
            }
        }
        if(code[7]) {
            try {
                BenchMark lj3 = new BenchMark();
                lj3.set(objectives, matrixA, b, lower, upper, nsize, directions, 4);
                ArrayList<Object> results = lj3.solve();
                main.useResults(results, nsize, "BenchMark4", output);
                xb = (ArrayList<Double>) results.get(0);
            } catch (IloException e) {
                e.printStackTrace();
            }
        }
        if(code[8]) {
            try {
                CLJ7 lj3 = new CLJ7();
                lj3.set(objectives, matrixA, b, lower, upper, nsize, directions);
                ArrayList<Object> results = lj3.solve();
                //main.useResults(results, nsize, "CLJ7          ");
                xc = (ArrayList<Double>) results.get(0);
            } catch (IloException e) {
                e.printStackTrace();
            }

        }
        if(code[9]) {
            try {
                CLJ8 lj3 = new CLJ8();
                lj3.set(objectives, matrixA, b, lower, upper, nsize);
                ArrayList<Object> results = lj3.solve();
                //main.useResults(results, nsize, "CLJ8          ");
            } catch (IloException e) {
                e.printStackTrace();
            }

        }
        if(code[10]) {
            try {
                CLJ9 lj3 = new CLJ9();
                lj3.set(objectives, matrixA, b, lower, upper, nsize, directions);
                ArrayList<Object> results = lj3.solve();
                //main.useResults(results, nsize, "CLJ8          ");
            } catch (IloException e) {
                e.printStackTrace();
            }

        }

        //System.out.println(" ////////////// chebyDiff " + FButils.L1norm(xb,xc,intvar));
    }

    public static ArrayList<Object> RANDOMINSIDEABOX( int intvar, int convar, int pp){
        int m = 1;
        int n = intvar + convar;//nsize;// (int) Math.round(Math.random()*200);
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
            b[i] =(10 - Math.random()*5)*100;
            //System.out.println(b[i]);
        }

        double valup = 100;
        double vallow = 0;


        for(int i =0 ; i< intvar; i++){
            lower[i] = vallow;
            upper[i] = valup;
        }
        for(int i = intvar; i < n; i++){
            lower[i] = 0;
            upper[i] = 10;
        }
        //////////////////////////////////////




        for(int i = 0; i<m; i++){
            for(int j = 0; j<n; j++)
                o[i][j] = Math.round(Math.random()*10000);//-50;
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

    private static void useResults(ArrayList<Object> results, int nsize, String algorithm_name, String output) {
        ArrayList<Double> x = (ArrayList<Double>) results.get(0);
        double [] ret = new double [6];
        double ig = FButils.integralityGap(IntegralityGapTypes.L2Norm, nsize, x);
        ret[0] = ig;
        boolean stop = (boolean) results.get(1);
        ret[1] = (stop? 0d : 1d);
        int iter = (int) results.get(2);
        ret[2] = iter;
        long time = (long) results.get(3);
        ret[3] = time;
        boolean check = (boolean) results.get(4);
        ret[4] = (check? 0d : 1d);
        double objgap = (double) results.get(5);
        ret[5] = objgap;
        System.out.println(algorithm_name+"  IntegralityGap: " + ig + "  Stop: " + stop + "  Iter: " + iter + "  Time (msec): " + (time ) + "   Check:" + check + "  ObjGap: " + objgap);


        try {
            boolean verify = false;
            if (!new File(output).exists()) {
                verify = true;
            }
            FileWriter file = new FileWriter(output, !verify);
            CsvWriter outputWriterkh = new CsvWriter(file, ',');
            if (verify) {
                outputWriterkh.write("Name");

                outputWriterkh.write("Integrality_Gap");
                outputWriterkh.write("Stop_Condition");
                outputWriterkh.write("Iterations");
                outputWriterkh.write("Time(msec)");
                outputWriterkh.write("Check_Feasibility");
                outputWriterkh.write("Objvalue");
                outputWriterkh.endRecord();
            }

            outputWriterkh.write(algorithm_name);
            for (int j = 0; j < ret.length; j++) {
                outputWriterkh.write(ret[j] + "");
            }
            outputWriterkh.endRecord();
            outputWriterkh.close();
        }catch(IOException e){
            e.printStackTrace();
        }

    }



    private static void printMatrix(double [] p){
        System.out.println(".............................");
        for(int i = 0 ; i < p.length; i++)
            System.out.println( i+" >  "+p[i]);
    }
    private static void printMatrix(int [] p){
        System.out.println(".............................");
        for(int i = 0 ; i < p.length; i++)
            System.out.println( i+" >  "+p[i]);
    }
    private static void printMatrix(double [][] p){
        System.out.println(".............................");
        for(int i = 0 ; i < p.length; i++) {
            System.out.println();
            for (int j = 0; j < p[0].length; j++ ){
                System.out.print( "  " + p[i][j]);}
        }
        System.out.println();
    }
}
