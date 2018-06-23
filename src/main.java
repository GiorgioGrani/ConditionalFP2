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
    private static int i ;
    public static void main(String [] aaa){

try {
    run(aaa[0], aaa[1], aaa[2]);
}catch( OutOfMemoryError e ){
   e.printStackTrace();
}
    }
    private static int find (ArrayList<String> as, String s, String base){

//    if (as.size() > 0) {
//        int n = as.size();
//        String ref = as.get(n - 1);
//        int i = n - 1;
//        int low = 0;
//        int up = n - 1;
//
//        double l = Math.log(n) / Math.log(2);
//        int count = 0;
//
//        while (count < l) {
//            count++;
//            int comp = ref.compareTo(s);
//            if (comp == 0) return i;
//            if (comp < 0) {
//                int piv = i;
//                i = (int) Math.max(Math.round((i - low) / 2d - 0.1), 0);
//                up = piv;
//            }
//            if (comp > 0) {
//                int piv = i;
//                i = (int) Math.max(Math.round((up - i) / 2d + 0.1), 0);
//                low = piv;
//            }
//        }
//    }
//
//    as.add(base + main.i);
//    main.i = main.i + 1;
//    return as.size() -1;


        int i = 0;
        for(String si : as){
            if(si.equalsIgnoreCase(s)) return i;

            i++;
        }

        as.add(base + main.i);
        main.i = main.i +1;
        return as.size() - 1;
    }

    private static int find (ArrayList<String> as, String s){

        int i = 0;
        for(String si : as){
            if(si.equalsIgnoreCase(s)) return i;

            i++;
        }


        return  - 1;
    }



    public static ArrayList<Object> readLP(String lp) throws FileNotFoundException{
        ArrayList<String> binaries = new ArrayList<>();
        ArrayList<String> generals = new ArrayList<>();
        ArrayList<String> continuous = new ArrayList<>();
        ArrayList<String> constraints = new ArrayList<>();
        TreeMap<String,  TreeMap <String, Double>> A = new TreeMap<>();
        TreeMap<String,  TreeMap <String,TreeMap<String,Double>>> Aq = new TreeMap<>();
        TreeMap<String, Double> b = new TreeMap<>();
        TreeMap<String, TreeMap<String, Double>> cq = new TreeMap<>();
        TreeMap<String, Double> c = new TreeMap<>();
        TreeMap<String, Double> lower = new TreeMap<>();
        TreeMap<String, Double> upper = new TreeMap<>();
        TreeMap<String, Integer> directions = new TreeMap<>();
        String direction = "";
        int n = 0, ncont = 0, nbin = 0 , nint = 0;
        int eq = 0, geq = 0, leq = 0;
        main.i = 2;
        Scanner scan = new Scanner(new File(lp));
        scan.nextLine();
        scan.nextLine();
        scan.next();
        scan.nextInt();
        eq = scan.nextInt();
        geq = scan.nextInt();
        leq = scan.nextInt();
        scan.nextLine();
        scan.nextLine();
        scan.nextLine();
        scan.nextLine();
        scan.nextLine();
        scan.next();
        n = scan.nextInt();
        ncont = scan.nextInt();
        nbin = scan.nextInt();
        nint = scan.nextInt();
        scan.nextLine();
        scan.nextLine();
        scan.nextLine();
        scan.nextLine();
        scan.nextLine();
        scan.nextLine();
        direction = scan.next();
        System.out.println(direction);
        scan.next();
        String preslider = "1d";
        String slider = scan.next();
        String compare = scan.next();
        double val = 0d;

        char[] sel0 = compare.toCharArray();
        boolean flag = false;
        if(slider.equals("-")){
            if(sel0[0] == 'x' || sel0[0] == 'b' || sel0[0] == 'i'){
                val = -1d;
                slider = compare;
                //compare = scan.next();
            }else{
                val = -1d*(Double.parseDouble(compare));
                slider = scan.next();
                //compare= scan.next();
            }
        }else{
            if(sel0[0] == 'x' || sel0[0] == 'b' || sel0[0] == 'i'){
                val = 1d*Double.parseDouble(slider);
                slider = compare;
                //compare = scan.next();
            }else{
                val = 1d;
                flag = true;
            }
        }


        while( !slider.equalsIgnoreCase("Subject")){
            if(flag) {
                flag = false;
            }else{
                compare = scan.next();
            }
            if(compare.equals("+") || compare.equals("-") || compare.equals("Subject")){
                if(!preslider.equals("*")){
                    char[] sel = slider.toCharArray();
                    if(sel[sel.length-2] == '^'){
                        String name = "";
                        for(int i = 0; i < sel.length-2; i ++ ) name += sel[i];
                        slider = name;
                        if (sel[0] == 'x') {
                            if(cq.get(continuous.get(find(continuous, slider, "x"))) == null){
                                TreeMap<String, Double> sub = new TreeMap<>();
                                cq.put(continuous.get(find(continuous, slider, "x")), sub);
                            }
                            cq.get(continuous.get(find(continuous, slider, "x"))).put(continuous.get(find(continuous, slider, "x")), val);
                        }
                        if (sel[0] == 'b'){
                            if(cq.get(binaries.get(find(binaries, slider, "x"))) == null){
                                TreeMap<String, Double> sub = new TreeMap<>();
                                cq.put(binaries.get(find(binaries, slider, "x")), sub);
                            }
                            cq.get(binaries.get(find(binaries, slider, "x"))).put(binaries.get(find(binaries, slider, "x")), val);
                        }
                        if (sel[0] == 'i'){
                            if(cq.get(generals.get(find(generals, slider, "x"))) == null){
                                TreeMap<String, Double> sub = new TreeMap<>();
                                cq.put(generals.get(find(generals, slider, "x")), sub);
                            }
                            cq.get(generals.get(find(generals, slider, "x"))).put(generals.get(find(generals, slider, "x")), val);
                        }
                    }else {
                        if (sel[0] == 'x') c.put(continuous.get(find(continuous, slider, "x")), val);
                        if (sel[0] == 'b') c.put(binaries.get(find(binaries, slider, "b")), val);
                        if (sel[0] == 'i') c.put(generals.get(find(generals, slider, "i")), val);
                    }
                }else{
                    preslider ="";
                }
                if(compare.equals("+")){
                    compare = scan.next();
                    sel0 = compare.toCharArray();
                    if(sel0[0] == 'x' || sel0[0] == 'b' || sel0[0] == 'i'){
                        val = 1d;
                        slider = compare;
                    }else{
                        val = Double.parseDouble(compare);
                        compare = scan.next();
                        slider = compare;
                    }
                }else if(compare.equals("-")) {
                    compare = scan.next();
                    sel0 = compare.toCharArray();
                    if(sel0[0] == 'x' || sel0[0] == 'b' || sel0[0] == 'i'){
                        val = -1d;
                        slider = compare;
                    }else{
                        val = -Double.parseDouble(compare);
                        compare = scan.next();
                        slider = compare;
                    }
                }
                else if(compare.equals("Subject"))  break;

            }else  if(compare.equals("*")){
                preslider = compare;
                String support = slider;
                compare = scan.next();
                sel0 = slider.toCharArray();
                char [] sel  = compare.toCharArray();
                //System.out.println(sel[0]+" "+sel0[0]);

                if(sel0[0] == 'x' && sel[0] == 'x') {
                    int i = find( continuous, slider, "x");
                    int j = find( continuous, compare, "x");

                    if(cq.get(continuous.get(i)) == null){
                        TreeMap<String, Double> sub = new TreeMap<>();
                        cq.put(continuous.get(i), sub);
                    }

                    cq.get(continuous.get(i)).put(continuous.get(j), val);
                }else if(sel0[0] == 'x' && sel[0] == 'b') {
                    int i = find( continuous, slider, "x");
                    int j = find( binaries, compare, "b");

                    if(cq.get(continuous.get(i)) == null){
                        TreeMap<String, Double> sub = new TreeMap<>();
                        cq.put(continuous.get(i), sub);
                    }

                    cq.get(continuous.get(i)).put(binaries.get(j), val);
                }else if(sel0[0] == 'b' && sel[0] == 'x') {
                    int i = find( binaries, slider, "b");
                    int j = find( continuous, compare, "x");

                    if(cq.get(binaries.get(i)) == null){
                        TreeMap<String, Double> sub = new TreeMap<>();
                        cq.put(binaries.get(i), sub);
                    }

                    cq.get(binaries.get(i)).put(continuous.get(j), val);
                }else if(sel0[0] == 'b' && sel[0] == 'i') {
                    int i = find( binaries, slider, "b");
                    int j = find( generals, compare, "i");

                    if(cq.get(binaries.get(i)) == null){
                        TreeMap<String, Double> sub = new TreeMap<>();
                        cq.put(binaries.get(i), sub);
                    }

                    cq.get(binaries.get(i)).put(generals.get(j), val);
                }else if(sel0[0] == 'i' && sel[0] == 'x') {
                    int i = find( generals, slider, "i");
                    int j = find( continuous, compare, "x");

                    if(cq.get(generals.get(i)) == null){
                        TreeMap<String, Double> sub = new TreeMap<>();
                        cq.put(generals.get(i), sub);
                    }

                    cq.get(generals.get(i)).put(continuous.get(j), val);
                }else if(sel0[0] == 'x' && sel[0] == 'i') {

                    int i = find( continuous, slider, "x");
                    int j = find( generals, compare, "i");


                    if(cq.get(continuous.get(i)) == null){
                        TreeMap<String, Double> sub = new TreeMap<>();
                        cq.put(continuous.get(i), sub);
                    }

                    cq.get(continuous.get(i)).put(generals.get(j), val);
                }else if(sel0[0] == 'b' && sel[0] == 'b') {
                    int i = find( binaries, slider, "b");
                    int j = find( binaries, compare, "b");


                    if(cq.get(binaries.get(i)) == null){
                        TreeMap<String, Double> sub = new TreeMap<>();
                        cq.put(binaries.get(i), sub);
                    }

                    cq.get(binaries.get(i)).put(binaries.get(j), val);
                }else if(sel0[0] == 'i' && sel[0] == 'b') {
                    int i = find( generals, slider, "i");
                    int j = find( binaries, compare,"b");



                    if(cq.get(generals.get(i)) == null){
                        TreeMap<String, Double> sub = new TreeMap<>();
                        cq.put(generals.get(i), sub);
                    }

                    cq.get(generals.get(i)).put(binaries.get(j), val);
                }else if(sel0[0] == 'i' && sel[0] == 'i') {
                    int i = find( generals, slider, "i");
                    int j = find( generals, compare , "i");


                    if(cq.get(generals.get(i)) == null){
                        TreeMap<String, Double> sub = new TreeMap<>();
                        cq.put(generals.get(i), sub);
                    }

                    cq.get(generals.get(i)).put(generals.get(j), val);
                }
            }

            // preslider = slider;
            //slider = scan.next();
        }
        scan.next(); //removing the word 'to'

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////


      String name = scan.next();
      while (!name.equals("Bounds") && !name.equals("Binary")  && !name.equals("General")) {
          //System.out.println(name);
          constraints.add(name);
          A.put(name, new TreeMap<String, Double>());
          Aq.put(name, new TreeMap<>());
          slider = scan.next();
          compare = scan.next();
          val = 0d;

          sel0 = compare.toCharArray();
          flag = false;
          if (slider.equals("-")) {
              if (sel0[0] == 'x' || sel0[0] == 'b' || sel0[0] == 'i') {
                  val = -1d;
                  slider = compare;
                  //compare = scan.next();
              } else {
                  val = -1d * (Double.parseDouble(compare));
                  slider = scan.next();
                  //compare= scan.next();
              }
          } else {
              if (sel0[0] == 'x' || sel0[0] == 'b' || sel0[0] == 'i') {
                  val = 1d * Double.parseDouble(slider);
                  slider = compare;
                  //compare = scan.next();
              } else {
                  val = 1d;
                  flag = true;
              }
          }
          //System.out.println(slider+flag);

          loop:
          while (fitOp(slider)) {
              //System.out.println(slider+flag);
              if (flag) {
                  flag = false;
              } else {

                 compare = scan.next();
              }
              //System.out.println("   "+compare+" "+slider+rfitOp(compare));
              if (compare.equals("+") || compare.equals("-") || rfitOp(compare)) {
                  if (!preslider.equals("*")) {
                      char[] sel = slider.toCharArray();
                      if(sel[sel.length-2] == '^'){
                          String namer = "";
                          for(int i = 0; i < sel.length-2; i ++ ) namer += sel[i];
                          slider = namer;
                          if (sel[0] == 'x') {
                              if(Aq.get(name).get(continuous.get(find(continuous, slider, "x"))) == null){
                                  TreeMap<String, Double> sub = new TreeMap<>();
                                  Aq.get(name).put(continuous.get(find(continuous, slider, "x")), sub);
                              }
                              Aq.get(name).get(continuous.get(find(continuous, slider, "x"))).put(continuous.get(find(continuous, slider, "x")), val);
                          }
                          if (sel[0] == 'b'){
                              if(Aq.get(name).get(binaries.get(find(binaries, slider, "x"))) == null){
                                  TreeMap<String, Double> sub = new TreeMap<>();
                                  Aq.get(name).put(binaries.get(find(binaries, slider, "x")), sub);
                              }
                              Aq.get(name).get(binaries.get(find(binaries, slider, "x"))).put(binaries.get(find(binaries, slider, "x")), val);
                          }
                          if (sel[0] == 'i'){
                              if(Aq.get(name).get(generals.get(find(generals, slider, "x"))) == null){
                                  TreeMap<String, Double> sub = new TreeMap<>();
                                  Aq.get(name).put(generals.get(find(generals, slider, "x")), sub);
                              }
                              Aq.get(name).get(generals.get(find(generals, slider, "x"))).put(generals.get(find(generals, slider, "x")), val);
                          }
                      }else {
                          if (sel[0] == 'x') A.get(name).put(continuous.get(find(continuous, slider, "x")), val);
                          if (sel[0] == 'b') A.get(name).put(binaries.get(find(binaries, slider, "b")), val);
                          if (sel[0] == 'i') A.get(name).put(generals.get(find(generals, slider, "i")), val);
                      }
                  } else {
                      preslider = "";
                  }
                  if (compare.equals("+")) {
                      compare = scan.next();
                      sel0 = compare.toCharArray();
                      if (sel0[0] == 'x' || sel0[0] == 'b' || sel0[0] == 'i') {
                          val = 1d;
                          slider = compare;
                      } else {
                          val = Double.parseDouble(compare);
                          compare = scan.next();
                          slider = compare;
                      }
                  } else if (compare.equals("-")) {
                      compare = scan.next();
                      sel0 = compare.toCharArray();
                      if (sel0[0] == 'x' || sel0[0] == 'b' || sel0[0] == 'i') {
                          val = -1d;
                          slider = compare;
                      } else {
                          val = -Double.parseDouble(compare);
                          compare = scan.next();
                          slider = compare;
                      }
                  } else if (rfitOp(compare)){ //todo primo cambio prima cera slider invc di compre

                      if(compare.equals("<=")){
                          directions.put(name, -1);
                      }else if (compare.equals(">=")){
                          directions.put(name, 1);
                      }else if (compare.equals("=")){
                          directions.put( name, 0);
                      }

                      String see = scan.next();
                      double track = Double.parseDouble(see);
                      b.put(name, track);


                      //System.out.println(compare+track+"  "+see);
                      name = scan.next();

                      break loop;
                  }

              } else if (compare.equals("*")) {
                  preslider = compare;
                  String support = slider;
                  compare = scan.next();
                  sel0 = slider.toCharArray();
                  char[] sel = compare.toCharArray();
                  if (sel0[0] == 'x' || sel[0] == 'x') {
                      int i = find(continuous, slider,"x");
                      int j = find(continuous, compare,"x");


                      if (Aq.get(name).get(continuous.get(i)) == null) {
                          TreeMap<String, Double> sub = new TreeMap<>();
                          Aq.get(name).put(continuous.get(i), sub);
                      }

                      Aq.get(name).get(continuous.get(i)).put(continuous.get(j), val);
                  } else if (sel0[0] == 'x' && sel[0] == 'b') {
                      int i = find(continuous, slider,"x");
                      int j = find(binaries, compare,"b");


                      if (Aq.get(name).get(continuous.get(i)) == null) {
                          TreeMap<String, Double> sub = new TreeMap<>();
                          Aq.get(name).put(continuous.get(i), sub);
                      }

                      Aq.get(name).get(continuous.get(i)).put(binaries.get(j), val);
                  }else if(sel0[0] == 'b' && sel[0] == 'x') {
                      int i = find( binaries, slider, "b");
                      int j = find( continuous, compare, "x");

                      if(Aq.get(name).get(binaries.get(i)) == null){
                          TreeMap<String, Double> sub = new TreeMap<>();
                          Aq.get(name).put(binaries.get(i), sub);
                      }

                      Aq.get(name).get(binaries.get(i)).put(continuous.get(j), val);
                  }else if(sel0[0] == 'b' && sel[0] == 'i') {
                      int i = find( binaries, slider, "b");
                      int j = find( generals, compare, "i");

                      if(Aq.get(name).get(binaries.get(i)) == null){
                          TreeMap<String, Double> sub = new TreeMap<>();
                          Aq.get(name).put(binaries.get(i), sub);
                      }

                      Aq.get(name).get(binaries.get(i)).put(generals.get(j), val);
                  }else if(sel0[0] == 'i' && sel[0] == 'x') {
                      int i = find( generals, slider, "i");
                      int j = find( continuous, compare, "x");

                      if(Aq.get(name).get(generals.get(i)) == null){
                          TreeMap<String, Double> sub = new TreeMap<>();
                          Aq.get(name).put(generals.get(i), sub);
                      }

                      Aq.get(name).get(generals.get(i)).put(continuous.get(j), val);
                  } else if (sel0[0] == 'x' && sel[0] == 'i') {
                      int i = find(continuous, slider,"x");
                      int j = find(generals, compare,"i");

                      if (Aq.get(name).get(continuous.get(i)) == null) {
                          TreeMap<String, Double> sub = new TreeMap<>();
                          Aq.get(name).put(continuous.get(i), sub);
                      }

                      Aq.get(name).get(continuous.get(i)).put(generals.get(j), val);
                  } else if (sel0[0] == 'b' && sel[0] == 'b') {
                      int i = find(binaries, slider,"b");
                      int j = find(binaries, compare,"b");

                      if (Aq.get(name).get(binaries.get(i)) == null) {
                          TreeMap<String, Double> sub = new TreeMap<>();
                          Aq.get(name).put(binaries.get(i), sub);
                      }

                      Aq.get(name).get(binaries.get(i)).put(binaries.get(j), val);
                  } else if (sel0[0] == 'i' && sel[0] == 'b') {
                      int i = find(generals, slider,"i");
                      int j = find(binaries, compare,"b");

                      if (Aq.get(name).get(generals.get(i)) == null) {
                          TreeMap<String, Double> sub = new TreeMap<>();
                          Aq.get(name).put(generals.get(i), sub);
                      }

                      Aq.get(name).get(generals.get(i)).put(binaries.get(j), val);
                  } else if (sel0[0] == 'i' && sel[0] == 'i') {
                      int i = find(generals, slider,"i");
                      int j = find(generals, compare,"i");

                      if (Aq.get(name).get(generals.get(i)) == null) {
                          TreeMap<String, Double> sub = new TreeMap<>();
                          Aq.get(name).put(generals.get(i), sub);
                      }

                      Aq.get(name).get(generals.get(i)).put(generals.get(j), val);
                  }
              }

          }

      }

      if(name.equals("Bounds")){
          slider = scan.next();
          while(!slider.equals("Binary") && !slider.equals("General") &&
                  !slider.equals("End") && !slider.equals("Ranges")){
              char[] sel = slider.toCharArray();
              compare = scan.next();
              val = Double.parseDouble(scan.next());
              if(sel[0] == 'x'){
                  if(compare.equals("<=")){
                      upper.put(continuous.get(find(continuous, slider,"x")), val);
                  }else if(compare.equals(">=")){
                      lower.put(continuous.get(find(continuous, slider,"x")), val);
                  }else if(compare.equals("=")){
                      upper.put(continuous.get(find(continuous, slider,"x")), val);
                      lower.put(continuous.get(find(continuous, slider,"x")), val);
                  }
              }else if (sel[0] == 'b'){
                  if(compare.equals("<=")){
                      upper.put(binaries.get(find(binaries, slider,"b")), val);
                  }else if(compare.equals(">=")){
                      lower.put(binaries.get(find(binaries, slider,"b")), val);
                  }else if(compare.equals("=")){
                      upper.put(binaries.get(find(binaries, slider,"b")), val);
                      lower.put(binaries.get(find(binaries, slider,"b")), val);
                  }
              }else if (sel[0] == 'i'){
                  if(compare.equals("<=")){
                      upper.put(generals.get(find(generals, slider,"i")), val);
                  }else if(compare.equals(">=")){
                      lower.put(generals.get(find(generals, slider,"i")), val);
                  }else if(compare.equals("=")){
                      upper.put(generals.get(find(generals, slider,"i")), val);
                      lower.put(generals.get(find(generals, slider,"i")), val);
                  }
              }
              slider = scan.next();
          }
      }



        ArrayList<Object> ret = new ArrayList<>();

        ret.add(n);
        ret.add(ncont);
        ret.add(nbin);
        ret.add(nint);
        ret.add(eq);
        ret.add(geq);
        ret.add(leq);
        ret.add(continuous);
        ret.add(binaries);
        ret.add(generals);
        ret.add(constraints);
        ret.add(direction);
        ret.add(c);
        ret.add(cq);
        ret.add(A);
        ret.add(Aq);
        ret.add(b);
        ret.add(directions);
        ret.add(lower);
        ret.add(upper);
        return ret;
    }

    private static boolean fitOp(String slider){
        return !slider.equalsIgnoreCase(">") && !slider.equalsIgnoreCase("<") && !slider.equalsIgnoreCase("<=") && !slider.equalsIgnoreCase(">=")
                && !slider.equalsIgnoreCase("=");
    }
    private static boolean rfitOp(String slider){
        return slider.equalsIgnoreCase(">") || slider.equalsIgnoreCase("<") || slider.equalsIgnoreCase("<=") || slider.equalsIgnoreCase(">=")
                || slider.equalsIgnoreCase("=");
    }

    public static ArrayList<Object> readMPS(String mps) throws FileNotFoundException{
//        ArrayList<String> binaries = new ArrayList<>();
//        ArrayList<String> generals = new ArrayList<>();
//        ArrayList<String> continuous = new ArrayList<>();
//        ArrayList<String> constraints = new ArrayList<>();
//        TreeMap<String,  TreeMap <String, Double>> A = new TreeMap<>();
//        TreeMap<String,  TreeMap <String,TreeMap<String,Double>>> Aq = new TreeMap<>();
//        TreeMap<String, Double> b = new TreeMap<>();
//        TreeMap<String, TreeMap<String, Double>> cq = new TreeMap<>();
//        TreeMap<String, Double> c = new TreeMap<>();
//        TreeMap<String, Double> lower = new TreeMap<>();
//        TreeMap<String, Double> upper = new TreeMap<>();
//        TreeMap<String, Integer> directions = new TreeMap<>();
//        String direction = "";
//        int n = 0, ncont = 0, nbin = 0 , nint = 0;
//        int eq = 0, geq = 0, leq = 0;

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
            while (!next.equalsIgnoreCase("BOUNDS") && !next.equalsIgnoreCase("ENDATA") &&  !next.equalsIgnoreCase("RANGES")) {
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

//        ArrayList<String> binaries = new ArrayList<>();
//        ArrayList<String> generals = new ArrayList<>();
//        ArrayList<String> continuous = new ArrayList<>();
//        ArrayList<String> constraints = new ArrayList<>();
//        TreeMap<String,  TreeMap <String, Double>> A = new TreeMap<>();
//        TreeMap<String,  TreeMap <String,TreeMap<String,Double>>> Aq = new TreeMap<>();
//        TreeMap<String, Double> b = new TreeMap<>();
//        TreeMap<String, TreeMap<String, Double>> cq = new TreeMap<>();
//        TreeMap<String, Double> c = new TreeMap<>();
//        TreeMap<String, Double> lower = new TreeMap<>();
//        TreeMap<String, Double> upper = new TreeMap<>();
//        TreeMap<String, Integer> directions = new TreeMap<>();
//        String direction = "";
//        int n = 0, ncont = 0, nbin = 0 , nint = 0;
//        int eq = 0, geq = 0, leq = 0;

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
    public static ArrayList<Object> readMPS2(String mps) throws FileNotFoundException{
        ArrayList<String> binaries = new ArrayList<>();
        ArrayList<String> continuous = new ArrayList<>();
        ArrayList<String> constraints = new ArrayList<>();
        TreeMap<String,  TreeMap <String, Double>> A = new TreeMap<>();
        TreeMap<String, Double> b = new TreeMap<>();
        TreeMap<String, Double> c = new TreeMap<>();
        TreeMap<String, Double> lower = new TreeMap<>();
        TreeMap<String, Double> upper = new TreeMap<>();
        TreeMap<String, Integer> directions = new TreeMap<>();
        String direction = "";
        int n = 0, ncont = 0, nbin = 0 , nint = 0;
        int eq = 0, geq = 0, leq = 0;

        Scanner scan = new Scanner(new File(mps));
        String name = mps;
        String objname = "obj";

        String next = scan.next();
        if(next.equalsIgnoreCase("NAME")){
            name = scan.next();
            next = scan.next();
        }
        if(next.equalsIgnoreCase("ROWS")){
            next = scan.next();
            while( !next.equalsIgnoreCase("COLUMNS")){
                String directiond = next;
                String conname = scan.next();
if(directiond.equalsIgnoreCase("N")) objname = conname;
                constraints.add(conname);
                if(directiond.equalsIgnoreCase("L")) directions.put(constraints.get(find(constraints,conname)), -1);
                if(directiond.equalsIgnoreCase("G")) directions.put(constraints.get(find(constraints,conname)),  1);
                if(directiond.equalsIgnoreCase("E")) directions.put(constraints.get(find(constraints,conname)),  0);
                next = scan.next();
            }
        }

        for(String con : constraints){
            if(!con.equalsIgnoreCase(objname)){
                A.put(con, new TreeMap<String, Double>());
            }
        }

        if (next.equalsIgnoreCase("COLUMNS")) {
            next = scan.next();
            while (!next.equalsIgnoreCase("RHS") && !next.equalsIgnoreCase("BOUNDS")) {
                if (next.equalsIgnoreCase("MARK0000") || next.equalsIgnoreCase("MARK") ) {
                    scan.next();
                    next = scan.next();
                    next = scan.next();
                    while (!next.equalsIgnoreCase("MARK0001") &&
                            !next.equalsIgnoreCase("MARKEND")) {


                        String var = next;
                        binaries.add(var);
                        TreeMap<String, Double> submap = new TreeMap<>();

                        while (true) {
                            next = scan.next();
                            //System.out.println(next);
                            if (scan.hasNextDouble()) {
                                Double d = scan.nextDouble();
                                submap.put(next, d);
                                if(next.equalsIgnoreCase(objname)) c.put(var, d);
                                else A.get(constraints.get(find(constraints, next))).put(var, d);
                            } else if (next.equalsIgnoreCase(var)) {

                            } else {
                                break;
                            }
                        }


                    }

                    scan.next();
                    scan.next();
                    next = scan.next();
                }
                if (!next.equalsIgnoreCase("RHS")) {
                    while (!next.equalsIgnoreCase("RHS") && !next.equalsIgnoreCase("BOUNDS") && !next.equalsIgnoreCase("MARK0000") && !next.equalsIgnoreCase("MARK")) {


                        String var = next;
                        continuous.add(var);
                        TreeMap<String, Double> submap = new TreeMap<>();

                        while (true) {
                            next = scan.next();
                            if (scan.hasNextDouble()) {
                                Double d = scan.nextDouble();
                                submap.put(next, d);
                                if(next.equalsIgnoreCase(objname)) c.put(var, d);
                                else A.get(constraints.get(find(constraints, next))).put(var, d);
                                //System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@           "+d);
                            } else if (next.equalsIgnoreCase(var)) {

                            } else {
                                break;
                            }
                        }


                    }

                }

            }

        }


        TreeMap<String, Double> rhs = new TreeMap<>();
        if (next.equalsIgnoreCase("RHS")) {

            next = scan.next();
            while (!next.equalsIgnoreCase("BOUNDS") && !next.equalsIgnoreCase("ENDATA") &&  !next.equalsIgnoreCase("RANGES")) {
                //System.out.println("
                next = scan.next();
                double df = scan.nextDouble();
                rhs.put(next, df);
                b.put(constraints.get(find(constraints,next)), df);

                next = scan.next();
                if(scan.hasNextDouble()){
                    double d = scan.nextDouble();
                    rhs.put(next, d);
                    //System.out.println(constraints.get(find(constraints,next))+" baboom");
                    b.put(constraints.get(find(constraints,next)), d);
                    next = scan.next();
                }
                //
            }

            for(String con : constraints){
                if(b.get(con) == null) b.put(constraints.get(find(constraints,con)), 0d);
            }
        }


        if (next.equalsIgnoreCase("BOUNDS")) {
            //next = scan.next();
            while (scan.hasNext()) {
                next = scan.next();
                //System.out.println(next);
                if (next.equalsIgnoreCase("UP") || next.equalsIgnoreCase("UI")) {
                    scan.next();
                    String xname = scan.next();
                    int ind = find(binaries, xname);
                    double d = scan.nextDouble();
                    if(ind <= -1) {
                        ind = find(continuous, xname);
                        upper.put(continuous.get(ind),d);
                    }else{
                        upper.put(binaries.get(ind),d);
                    }

                } else if (next.equalsIgnoreCase("LO") || next.equalsIgnoreCase("LI")) {

                    scan.next();
                    String xname = scan.next();
                    int ind = find(binaries, xname);
                    double d = scan.nextDouble();
                    if(ind <= -1) {
                        ind = find(continuous, xname);
                        lower.put(continuous.get(ind),d);
                    }else{
                        lower.put(binaries.get(ind),d);
                    }

                } else if (next.equalsIgnoreCase("FX")) {
                    scan.next();
                    String xname = scan.next();
                    Double xval = scan.nextDouble();
                    int ind = find(binaries, xname);
                    //double d = scan.nextDouble();
                    if(ind <= -1) {
                        ind = find(continuous, xname);
                        upper.put(continuous.get(ind),xval);
                        lower.put(continuous.get(ind),xval);
                    }else{
                        upper.put(binaries.get(ind),xval);
                        lower.put(binaries.get(ind),xval);
                    }

                } else if (next.equals("FR")) {
                    scan.next();
                    String xname = scan.next();

                    int ind = find(binaries, xname);
                    if(ind <= -1) {
                        ind = find(continuous, xname);
                        upper.put(continuous.get(ind),Double.MAX_VALUE);
                        lower.put(continuous.get(ind), -(Double.MAX_VALUE-1));
                    }else{
                        upper.put(binaries.get(ind),Double.MAX_VALUE);
                        lower.put(binaries.get(ind), -(Double.MAX_VALUE-1));
                    }

                } else if (next.equalsIgnoreCase("BV")) {
                    scan.next();
                    String xname = scan.next();
                    int ind = find(binaries, xname);
                    //double d = scan.nextDouble();
                    if(ind <= -1) {
                        ind = find(continuous, xname);
                        upper.put(continuous.get(ind),1d);
                        lower.put(continuous.get(ind),0d);
                    }else{
                        upper.put(binaries.get(ind),1d);
                        lower.put(binaries.get(ind),0d);
                    }

                } else if (next.equalsIgnoreCase("MI")) {
                    scan.next();
                    String xname = scan.next();
                    int ind = find(binaries, xname);
                    //double d = scan.nextDouble();
                    if(ind <= -1) {
                        ind = find(continuous, xname);
                        upper.put(continuous.get(ind),0d);
                        lower.put(continuous.get(ind),-(Double.MAX_VALUE-1));
                    }else{
                        upper.put(binaries.get(ind),0d);
                        lower.put(binaries.get(ind),0d);
                    }

                }else if (next.equalsIgnoreCase("PL")) {
                    scan.next();
                    String xname = scan.next();
                    int ind = find(binaries, xname);
                    //double d = scan.nextDouble();
                    if(ind <= -1) {
                        ind = find(continuous, xname);
                        upper.put(continuous.get(ind),Double.MAX_VALUE);
                        lower.put(continuous.get(ind),0d);
                    }else{
                        upper.put(binaries.get(ind),1d);
                        lower.put(binaries.get(ind),0d);
                    }
                    //System.out.println(xname+" yyyyyyyyyyyyyyyyyyyyyy  "+xval);
                }else if(next.equalsIgnoreCase("ENDATA")){
                    break;
                }
            }
        }


        constraints.remove(find(constraints,objname));
        n = binaries.size() + continuous.size();
        nint = 0 ;
        nbin = binaries.size();
        ncont = continuous.size();
        ArrayList<Object> ret = new ArrayList<>();

//        boolean ct = false, bt = false;
//        for (String s : binaries) if(s.equals("'INTEND'")) bt = true;
//        for (String s : continuous) if(s.equals("'INTEND'")) ct = true;
//        if(bt) binaries.remove(find(binaries, "'INTEND'"));
//        if(ct) continuous.remove(find(continuous, "'INTEND'"));
//        ct = false;
//        bt = false;
//        for (String s : binaries) if(s.equals("'INTORG'")) bt = true;
//        for (String s : continuous) if(s.equals("'INTORG'")) ct = true;
//        if(bt) binaries.remove(find(binaries, "'INTORG'"));
//        if(ct) continuous.remove(find(continuous, "'INTORG'"));
//        ct = false;
//        bt = false;
//        for (String s : binaries) if(s.equals("'INTORG'")) bt = true;
//        for (String s : continuous) if(s.equals("'INTORG'")) ct = true;
//        if(bt) binaries.remove(find(binaries, "'INTORG'"));
//        if(ct) continuous.remove(find(continuous, "'INTORG'"));

        ret.add(n);
        ret.add(ncont);
        ret.add(nbin);
        ret.add(nint);
        ret.add(eq);
        ret.add(geq);
        ret.add(leq);
        ret.add(continuous);
        ret.add(binaries);
        ret.add(new ArrayList<String>());
        ret.add(constraints);
        ret.add("minimize");
        ret.add(c);
        ret.add(new TreeMap<String, TreeMap<String, Double>>());
        ret.add(A);
        ret.add(new TreeMap<String, TreeMap<String, TreeMap<String,Double>>>());
        ret.add(b);
        ret.add(directions);
        ret.add(lower);
        ret.add(upper);
        return ret;
    }

    public static int codify(String s){
        if(s.equalsIgnoreCase("L")) return -1;
        else if(s.equalsIgnoreCase("E")) return 0;
        else return 1;
    }

    public static void run( String mps, String output, String name){

        ArrayList<Object> param = new ArrayList<>();
        char[] let = mps.toCharArray();
        int siz = let.length;
        String type = ""+let[siz -3]+let[siz - 2]+ let[siz - 1];
        try {
            if(type.equalsIgnoreCase("mps")) param = readMPS2(mps);
            else if (type.equalsIgnoreCase(".lp")) param = readLP(mps);
        }catch (FileNotFoundException e){
            e.printStackTrace();
            return;
        }

        System.out.println("matrices stored");


        if(type.equalsIgnoreCase("mps")) {
            boolean [] code = {true,         true,   true,      true,false,true,true,true,false , false, false};

            if(code[2]){
                DADUSHT fp = new DADUSHT();
                try {
                    fp.set(param);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "DADUSHT ", output, name);
                }catch(IloException e){
                    e.printStackTrace();
                }
            }

            if(code[0]){
                BenchMarkLP fp = new BenchMarkLP();
                try {
                    fp.set(param, 1);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "BENCHMARK  ", output, name);
                }catch(IloException e){
                    e.printStackTrace();
                }
                System.gc();
            }
            if(code[0]){
                BenchMarkLP fp = new BenchMarkLP();
                try {
                    fp.set(param, 2);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "BENCHMARK2 ", output, name);
                }catch(IloException e){
                    e.printStackTrace();
                }
                System.gc();
            }
            if(code[1]){
                BBFPT fp = new BBFPT();
                try {
                    fp.set(param);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "BBFPT ", output, name);
                }catch(IloException e){
                    e.printStackTrace();
                }
                System.gc();
            }

//            nsize = (int) param.get(0);
//
//            objectives = (double[]) param.get(1);
//            matrixA = (double[][]) param.get(2);
//            b = (double[]) param.get(3);
//            lower = (double[]) param.get(4);
//            upper = (double[]) param.get(5);
//            directions = (int[]) param.get(6);
//
//
////        printMatrix(objectives);
////        printMatrix(matrixA);
////        printMatrix(b);
////        printMatrix(lower);
////        printMatrix(upper);
////        printMatrix(directions);
//
//        boolean [] code = {true,false,false,true,false,true,true,true,false , false, false};
//        if(code[0]){
//            BasicVersion fp = new BasicVersion();
//            try {
//                fp.set(objectives, matrixA, b, lower, upper, nsize, directions);
//                ArrayList<Object> results = fp.solve();
//                main.useResults(results, nsize, "Simple", output, name);
//            }catch(IloException e){
//                e.printStackTrace();
//            }
//        }
//        if(code[1]) {
//            ReinforcementFP fp = new ReinforcementFP();
//            try {
//                fp.set(objectives, matrixA, b, lower, upper, nsize, directions);
//                ArrayList<Object> results = fp.solve();
//                main.useResults(results, nsize, "Reinforcement", output, name);
//            }catch(IloException e){
//                e.printStackTrace();
//            }
//        }
//        if(code[2]) {
//            AggressiveRFP fp = new AggressiveRFP();
//            try {
//                fp.set(objectives, matrixA, b, lower, upper, nsize, directions);
//                ArrayList<Object> results = fp.solve();
//                main.useResults(results, nsize, "Aggressive", output, name);
//            }catch(IloException e){
//                e.printStackTrace();
//            }
//        }
//
//
//        if(code[3]) {
//
//                DADUSH fp = new DADUSH();
//            try {
//                fp.set(objectives, matrixA, b, lower, upper, nsize, directions);
//                ArrayList<Object> results = fp.solve();
//                main.useResults(results, nsize, "DADUSH", output, name);
//            }catch(IloException e){
//                e.printStackTrace();
//            }
//        }
//        if(code[5]) {
//            try {
//                BenchMark lj3 = new BenchMark();
//                lj3.set(objectives, matrixA, b, lower, upper, nsize, directions, 2);
//                ArrayList<Object> results = lj3.solve();
//                main.useResults(results, nsize, "BenchMark2", output, name);
//            } catch (IloException e) {
//                e.printStackTrace();
//            }
//        }
//        if(code[6]) {
//            try {
//                BenchMark lj3 = new BenchMark();
//                lj3.set(objectives, matrixA, b, lower, upper, nsize, directions,3 );
//                ArrayList<Object> results = lj3.solve();
//                main.useResults(results, nsize, "BenchMark3", output, name);
//            } catch (IloException e) {
//                e.printStackTrace();
//            }
//        }
//        if(code[7]) {
//            try {
//                BenchMark lj3 = new BenchMark();
//                lj3.set(objectives, matrixA, b, lower, upper, nsize, directions, 4);
//                ArrayList<Object> results = lj3.solve();
//                main.useResults(results, nsize, "BenchMark4", output, name);
//            } catch (IloException e) {
//                e.printStackTrace();
//            }
//        }
//        if(code[4]) {
//
//            BBFP fp = new BBFP();
//            try {
//                fp.set(objectives, matrixA, b, lower, upper, nsize, directions);
//                ArrayList<Object> results = fp.solve();
//                main.useResults(results, nsize, "BBFP", output, name);
//            }catch(IloException e){
//                e.printStackTrace();
//            }
//        }
//        if(code[8]) {
//            try {
//                CLJ7 lj3 = new CLJ7();
//                lj3.set(objectives, matrixA, b, lower, upper, nsize, directions);
//                ArrayList<Object> results = lj3.solve();
//                //main.useResults(results, nsize, "CLJ7
//            } catch (IloException e) {
//                e.printStackTrace();
//            }
//
//        }
//        if(code[9]) {
//            try {
//                CLJ8 lj3 = new CLJ8();
//                lj3.set(objectives, matrixA, b, lower, upper, nsize);
//                ArrayList<Object> results = lj3.solve();
//                //main.useResults(results, nsize, "CLJ8          ");
//            } catch (IloException e) {
//                e.printStackTrace();
//            }
//
//        }
//        if(code[10]) {
//            try {
//                CLJ9 lj3 = new CLJ9();
//                lj3.set(objectives, matrixA, b, lower, upper, nsize, directions);
//                ArrayList<Object> results = lj3.solve();
//                //main.useResults(results, nsize, "CLJ8          ");
//            } catch (IloException e) {
//                e.printStackTrace();
//            }
//
//        }
//
//        //System.out.println(" ////////////// chebyDiff " + FButils.L1norm(xb,xc,intvar));
        }else if(type.equalsIgnoreCase(".lp")){
             //                 benchmark2   bbfpt   dadusht
            boolean [] code = {true,         false,   true,      true,false,true,true,true,false , false, false};
            if(code[0]){
                BenchMarkLP fp = new BenchMarkLP();
                try {
                    fp.set(param, 2);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "BENCHMARK ", output, name);
                }catch(IloException e){
                    e.printStackTrace();
                }
            }
            if(code[1]){
                BBFPT fp = new BBFPT();
                try {
                    fp.set(param);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "BBFPT ", output, name);
                }catch(IloException e){
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
            }if(code[2]){
                DADUSHT fp = new DADUSHT();
                try {
                    fp.set(param);
                    ArrayList<Object> results = fp.solve();
                    main.useResults(results, "DADUSHT ", output, name);
                }catch(IloException e){
                    e.printStackTrace();
                }
            }
        }
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

    public static void useResults(ArrayList<Object> results, int nsize, String algorithm_name, String output, String problemname) {
        ArrayList<Double> x = (ArrayList<Double>) results.get(0);
        double [] ret = new double [6];
        double ig = FButils.integralityGap(IntegralityGapTypes.L1Norm, nsize, x);
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
                outputWriterkh.write("Problem");
                outputWriterkh.write("Name");

                outputWriterkh.write("Integrality_Gap");
                outputWriterkh.write("Stop_Condition");
                outputWriterkh.write("Iterations");
                outputWriterkh.write("Time(msec)");
                outputWriterkh.write("Check_Feasibility");
                outputWriterkh.write("Objvalue");
                outputWriterkh.endRecord();
            }

            outputWriterkh.write(problemname);
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

    public static void useResults(ArrayList<Object> results, String algorithm_name, String output, String problemname) {
        TreeMap<String, Double> x = (TreeMap<String, Double>) results.get(0);
        double [] ret = new double [6];
        double ig = FButils.integralityGap(IntegralityGapTypes.L2Norm, x);
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
                outputWriterkh.write("Problem");
                outputWriterkh.write("Name");

                outputWriterkh.write("Integrality_Gap");
                outputWriterkh.write("Stop_Condition");
                outputWriterkh.write("Iterations");
                outputWriterkh.write("Time(msec)");
                outputWriterkh.write("Check_Feasibility");
                outputWriterkh.write("Objvalue");
                outputWriterkh.endRecord();
            }

            outputWriterkh.write(problemname);
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
