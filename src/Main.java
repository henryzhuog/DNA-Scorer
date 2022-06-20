import java.io.FileNotFoundException;
import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

public class Main {
    // BLOSUM62 Scoring Matrix
    public static final int[][] BLOSUM62 = {
        //   A   R   N   D   C   Q   E   G   H  I   L   K   M   F   P   S   T   W  Y  V
    /**A**/ {4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0},
    /**R**/ {-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3},
    /**N**/ {-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3},
    /**D**/ {-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3},
    /**C**/ {0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
    /**Q**/ {-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2},
    /**E**/ {-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2},
    /**G**/ {0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3},
    /**H**/ {-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3},
    /**I**/ {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3},
    /**L**/ {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1},
    /**K**/ {-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2},
    /**M**/ {-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1},
    /**F**/ {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1},
    /**P**/ {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2},
    /**S**/ {1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2},
    /**T**/ {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0},
    /**W**/ {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3},
    /**Y**/ {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1},
    /**V**/ {0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4}
    };

    // Order of nucleotides based on BLOSUM62 Scoring Matrix
    public static final char[] letterOrder = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H',
            'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
            
    // Names of FASTA files to be parsed
    public static final String[] testOrder = {"O95363.fasta.txt", "P10085.fasta.txt",
            "P13904.fasta.txt", "P15172.fasta.txt", "P16075.fasta.txt", "P17542.fasta.txt",
            "P22816.fasta.txt", "Q8IU24.fasta.txt", "Q10574.fasta.txt", "Q90477.fasta.txt"};

    // Gap penalty
    public static final int gapCost = -4;

    // Information pertaining to sequences 1 and 2
    public static Stack<Character> st1 = new Stack<Character>();
    public static Stack<Character> st2 = new Stack<Character>();
    public static String subSq1, subSq2;

    // x & y "correction"
    public static int xCorrection = 1;
    public static int yCorrection = 1;

    public static void main(String[] args) throws FileNotFoundException {
      // Uncomment lines 58, 59, 64, 67 119, & 120 
      // to run multiple iterations for increased accuracy
      
      //  for (int q = 0; q < 10; q++) {
      //      for (int w = 0; w < 10; w++) {
                Scanner s = new Scanner(System.in);
                // File name is "____.fasta.txt"
                System.out.print("Enter the file containing s1: ");
                String s1Name = s.nextLine();
                // String s1Name = testOrder[q];
                System.out.print("Enter the file containing s2: ");
                String s2Name = s.nextLine();
                // String s2Name = testOrder[w];
                System.out.print("Enter the number of permutations, N (must be > 0): ");
                int N = s.nextInt();
                System.out.println();

                File f1 = new File(s1Name);
                Scanner f1Reader = new Scanner(f1);
                File f2 = new File(s2Name);
                Scanner f2Reader = new Scanner(f2);

                String s1 = "";
                String s2 = "";
                String id1 = s1Name.replace(".fasta.txt", "");
                String id2 = s2Name.replace(".fasta.txt", "");

                s1 = parse(f1Reader, s1);
                s2 = parse(f2Reader, s2);

                int[][] SWMatrix = new int[s1.length() + 1][s2.length() + 1];

                buildSW(SWMatrix, s1, s2);

                int[] OPInfo = getOPInfo(SWMatrix, s1, s2);
                // 0-based, max int coordinate
                int iInt = OPInfo[1]; // y direction (sic)
                int jInt = OPInfo[2]; // x direction (sic)

                traceback(SWMatrix, s1, s2, iInt, jInt);

                // Form the aligning subsequence that is to be printed to output
                ArrayList<Character> arr1 = new ArrayList<>();
                ArrayList<Character> arr2 = new ArrayList<>();
                while (!st1.empty() || !st2.empty()) {
                    arr1.add(st1.pop());
                    arr2.add(st2.pop());
                }
                subSq1 = ALtoString(arr1);
                subSq2 = ALtoString(arr2);

                // Split into strings of length 60
                ArrayList<String> shortSQ1 = getShortSQ(subSq1);
                ArrayList<String> shortSQ2 = getShortSQ(subSq2);
                for (int i = 0; i < shortSQ1.size(); i++) {
                    System.out.println(id1 + ": " + (iInt + yCorrection + (i * 60)) +
                            " " + shortSQ1.get(i));
                    System.out.println();
                    System.out.println(id2 + ": " + (jInt + xCorrection + (i * 60)) +
                            " " + shortSQ2.get(i));
                    System.out.println();
                }
                System.out.println("Score of Optimal Alignment: " + OPInfo[0]);
                System.out.println();
        //    }
        // }

        // Print matrix if both s1 and s2 are shorter than 15 characters
        if (s1.length() < 15 && s2.length() < 15) {
            for (int[] row : SWMatrix)
                System.out.println(Arrays.toString(row));
        } else {
            System.out.println();
            System.out.println("Matrix is too big to print. One or more strings is longer" +
                    " than 15 characters.");
        }
        
        int[] permScores = new int[N];
        double kValue = 0.0;
        
        // Run the algorithm on a scrambled version of sequence 2
        for (int i = 0; i < permScores.length; i++) {
            String perm = rearrangeString(s2);
            int[][] permSW = new int[s1.length() + 1][perm.length() + 1];
            buildSW(permSW, s1, perm);
            int[] permInfo = getOPInfo(permSW, s1, perm);
            if (permInfo[0] >= OPInfo[0])
                kValue++;
            permScores[i] = permInfo[0];
        }

        double pValue = (kValue + 1.0)/(N + 1.0);
        System.out.println();
        NumberFormat numFormat = new DecimalFormat();
        numFormat = new DecimalFormat("0.######E0");
        System.out.println("p value = " + numFormat.format(pValue));
    }

    // Takes a DNA sequence s, scrambles it and returns the result
    public static String rearrangeString(String s) {
        char[] cArray = s.toCharArray();
        Random r = new Random();
        for (int i = cArray.length - 1; i > 0; i--) {
            int j = r.nextInt(i);
            char temp = cArray[i];
            cArray[i] = cArray[j];
            cArray[j] = temp;
        }
        String result = "";
        for (char c:cArray) {
            result += c;
        }
        return result;
    }

    // Takes an ArrayList of Characters arr, returns a String consisting
    // of all the characters in that ArrayList concatenated
    public static String ALtoString(ArrayList<Character> arr) {
        String s = "";
        for (Character a : arr) {
            s += a.toString();
        }
        return s;
    }

    // Takes a DNA sequence s, returns an ArrayList of Strings where each String is a 
    // 60 nucleotide k-mer of s
    public static ArrayList<String> getShortSQ(String s) {
        ArrayList<String> result = new ArrayList<String>();
        int index = 0;
        while (index < s.length()) {
            result.add(s.substring(index, Math.min(index + 60,s.length())));
            index += 60;
        }
        return result;
    }

    // Parse the sequences from the desired FASTA file
    public static String parse(Scanner sc, String s) {
        while (sc.hasNextLine()) {
            String data = sc.nextLine();
            if (data.charAt(0) != '>') {
                s += data;
            }
        }
        return s;
    }

    // Build the SW Matrix
    public static void buildSW(int[][] SWMatrix, String s1, String s2) {
        // Base cases, ref. slide 41
        SWMatrix[0][0] = 0;
        // First row
        for (int i = 1; i <= s1.length(); i++) {
            SWMatrix[i][0] = 0;
        }
        // First column
        for (int j = 1; j <= s2.length(); j++) {
            SWMatrix[0][j] = 0;
        }
        // Everything else
        for (int i = 1; i <= s1.length(); i++) {
            for (int j = 1; j <= s2.length(); j++) {
                /**
                 *      V(i,j) = max {  1. V(i-1,j-1) + match/mismatch  (diagonal)
                 *                      2. V(i-1, j) + gap (left to right, -->)
                 *                      3. V(i, j-1) + gap (down, \/)
                 *                      4. 0
                 *                   }
                 */
                char c1 = s1.charAt(i - 1);
                char c2 = s2.charAt(j - 1);
                // matchScore uses BLOSUM62 to determine the appropriate amount to add.
                int diag = SWMatrix[i - 1][j - 1] + matchScore(c1, c2);
                int horiz = SWMatrix[i-1][j] + gapCost; //O(S[i], ~) = gap
                int vert = SWMatrix[i][j-1] + gapCost; //O(~,T[j]) = gap
                SWMatrix[i][j] = Math.max(diag, Math.max(horiz, Math.max(vert, 0)));
            }
        }
    }

    // Uses BLOSUM62 to return the appropriate score
    public static int matchScore(char c1, char c2) {
        int a = 0;
        for (int i = 0; i < letterOrder.length; i++) {
            if (letterOrder[i] == c1) {
                a = i;
            }
        }
        int b = 0;
        for (int j = 0; j < letterOrder.length; j++) {
            if (letterOrder[j] == c2) {
                b = j;
            }
        }
        return BLOSUM62[a][b];
    }

    // Return an array containing the OP score and its matrix coordinate
    public static int[] getOPInfo(int[][] SWMatrix, String s1, String s2) {
        int[] result = new int[3];
        for (int i = 1; i <= s1.length(); i++) {
            for (int j = 1; j <= s2.length(); j++) {
                if (SWMatrix[i][j] > result[0]) {
                    result[0] = SWMatrix[i][j];
                    result[1] = i;
                    result[2] = j;
                }
            }
        }
        return result;
    }

    // Traceback function based on lecture slide
    public static void traceback(int[][] SWMatrix, String s1, String s2, int i, int j) {
        // exit on 0
        if (SWMatrix[i][j] == 0) {
            return;
        }
        // Traces back the highest int until it gets to 0
        if (SWMatrix[i][j] == SWMatrix[i-1][j] + gapCost) {
            yCorrection--;
            st1.add(s1.charAt(i-1));
            st2.add('-');
            traceback(SWMatrix, s1, s2, i - 1, j);
        } else if(SWMatrix[i][j] == SWMatrix[i][j-1] + gapCost) {
            xCorrection--;
            st1.add('-');
            st2.add(s2.charAt(j-1));
            traceback(SWMatrix, s1, s2, i, j - 1);
        } else {
            yCorrection--;
            xCorrection--;
            st1.push(s1.charAt(i - 1));
            st2.push(s2.charAt(j - 1));
            traceback(SWMatrix, s1, s2, i - 1, j - 1);
        }
    }
}