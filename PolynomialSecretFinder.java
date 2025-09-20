import java.io.*;
import java.math.BigInteger;
import java.util.*;
import com.google.gson.*;

public class PolynomialSecretFinder {

    static class Point {
        int x;
        BigInteger y;
        Point(int x, BigInteger y) {
            this.x = x;
            this.y = y;
        }
    }

    public static void main(String[] args) throws Exception {
        if(args.length < 1) {
            System.out.println("Usage: java -cp .;gson-2.10.1.jar PolynomialSecretFinder <input.json>");
            return;
        }

        String jsonText = new String(java.nio.file.Files.readAllBytes(java.nio.file.Paths.get(args[0])));
        JsonObject json = JsonParser.parseString(jsonText).getAsJsonObject();

        JsonObject keys = json.getAsJsonObject("keys");
        int n = keys.get("n").getAsInt();
        int k = keys.get("k").getAsInt();
        int degree = k - 1;

        List<Point> points = new ArrayList<>();
        for(Map.Entry<String, JsonElement> entry : json.entrySet()) {
            String key = entry.getKey();
            if(key.equals("keys")) continue;

            int x = Integer.parseInt(key);
            JsonObject valObj = entry.getValue().getAsJsonObject();
            int base = Integer.parseInt(valObj.get("base").getAsString());
            String value = valObj.get("value").getAsString();

            BigInteger y = decodeValue(value, base);
            points.add(new Point(x, y));
        }

        // Use first k points
        List<Point> chosen = points.subList(0, k);

        BigInteger secretC = solveMatrix(chosen, degree);
        System.out.println("Secret constant c = " + secretC);
    }

    // Decode value from any base
    static BigInteger decodeValue(String val, int base) {
        BigInteger res = BigInteger.ZERO;
        for (char ch : val.toLowerCase().toCharArray()) {
            int digit;
            if (ch >= '0' && ch <= '9') digit = ch - '0';
            else digit = ch - 'a' + 10;
            if (digit >= base) throw new RuntimeException("Invalid digit for base");
            res = res.multiply(BigInteger.valueOf(base)).add(BigInteger.valueOf(digit));
        }
        return res;
    }

    // Solve using Gaussian elimination with BigInteger
    static BigInteger solveMatrix(List<Point> pts, int degree) {
        int n = pts.size();
        BigInteger[][] mat = new BigInteger[n][degree + 1];
        BigInteger[] rhs = new BigInteger[n];

        // Build matrix
        for(int i=0;i<n;i++) {
            BigInteger xPow = BigInteger.ONE;
            for(int j=degree;j>=0;j--) {
                mat[i][j] = xPow;
                xPow = xPow.multiply(BigInteger.valueOf(pts.get(i).x));
            }
            rhs[i] = pts.get(i).y;
        }

        // Gaussian elimination
        for(int i=0;i<=degree;i++) {
            // Find pivot
            int pivot = i;
            for(int j=i+1;j<n;j++)
                if(mat[j][i].abs().compareTo(mat[pivot][i].abs()) > 0)
                    pivot = j;
            // Swap rows
            BigInteger[] tmp = mat[i]; mat[i]=mat[pivot]; mat[pivot]=tmp;
            BigInteger t = rhs[i]; rhs[i]=rhs[pivot]; rhs[pivot]=t;

            // Eliminate
            for(int j=i+1;j<n;j++) {
                BigInteger factor = mat[j][i].divide(mat[i][i]);
                for(int k=i;k<=degree;k++)
                    mat[j][k] = mat[j][k].subtract(factor.multiply(mat[i][k]));
                rhs[j] = rhs[j].subtract(factor.multiply(rhs[i]));
            }
        }

        // Back substitution
        BigInteger[] coeff = new BigInteger[degree+1];
        for(int i=degree;i>=0;i--) {
            BigInteger sum = rhs[i];
            for(int j=i+1;j<=degree;j++)
                sum = sum.subtract(mat[i][j].multiply(coeff[j]));
            coeff[i] = sum.divide(mat[i][i]);
        }

        return coeff[degree]; // constant term c
    }
}