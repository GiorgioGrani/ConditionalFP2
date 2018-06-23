import com.csvreader.CsvWriter;
import ilog.concert.IloException;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.TreeMap;
//import java.util.List;

public class mainCleaned {
    private static final double eps = 1e-8;
    private static int i;

    public static void main(String[] aaa) {

        try {
            run(aaa[0], aaa[1], aaa[2]);
        } catch (OutOfMemoryError e) {
            e.printStackTrace();
        }
    }

    public static void run(String mps, String output, String name) {

        ArrayList<Object> param = new ArrayList<>();
        char[] let = mps.toCharArray();
        int siz = let.length;
        String type = "" + let[siz - 3] + let[siz - 2] + let[siz - 1];
        try {
            if (type.equalsIgnoreCase("mps")) param = main.readMPS2(mps);
            else if (type.equalsIgnoreCase(".lp")) param = main.readLP(mps);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        }

        System.out.println("matrices stored");


        if (type.equalsIgnoreCase("mps")) {
            boolean[] code = {true, true, true, true, false, true, true, true, false, false, false};
            if (code[2]) {
                DADUSHT fp = new DADUSHT();
                try {
                    fp.set(param);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "DADUSHT ", output, name);
                } catch (IloException e) {
                    e.printStackTrace();
                }
            }

            if (code[0]) {
                BenchMarkLP fp = new BenchMarkLP();
                try {
                    fp.set(param, 1);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "BENCHMARK  ", output, name);
                } catch (IloException e) {
                    e.printStackTrace();
                }
                System.gc();
            }
            if (code[0]) {
                BenchMarkLP fp = new BenchMarkLP();
                try {
                    fp.set(param, 2);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "BENCHMARK2 ", output, name);
                } catch (IloException e) {
                    e.printStackTrace();
                }
                System.gc();
            }
            if (code[1]) {
                BBFPT fp = new BBFPT();
                try {
                    fp.set(param);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "BBFPT ", output, name);
                } catch (IloException e) {
                    e.printStackTrace();
                }
                System.gc();
            }


        } else if (type.equalsIgnoreCase(".lp")) {
            //                 benchmark2   bbfpt   dadusht
            boolean[] code = {true, false, true, true, false, true, true, true, false, false, false};
            if (code[0]) {
                BenchMarkLP fp = new BenchMarkLP();
                try {
                    fp.set(param, 2);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "BENCHMARK ", output, name);
                } catch (IloException e) {
                    e.printStackTrace();
                }
            }
            if (code[1]) {
                BBFPT fp = new BBFPT();
                try {
                    fp.set(param);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "BBFPT ", output, name);
                } catch (IloException e) {
                    e.printStackTrace();
                }
//            }if(code[0]){
//                BenchMarkLP fp = new BenchMarkLP();
//                try {
//                    fp.set(param, -1);
//                    ArrayList<Object> results = fp.solve();
//                    main.useResults(results, "BENCHMARK ", output, name);
//                }catch(IloException e){
//                    e.printStackTrace();
//                }
            }
            if (code[2]) {
                DADUSHT fp = new DADUSHT();
                try {
                    fp.set(param);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "DADUSHT ", output, name);
                } catch (IloException e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
