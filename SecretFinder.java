package org.example;// SecretFinder.java
import java.math.BigInteger;
import java.nio.file.*;
import java.util.*;
import com.google.gson.*;

public class SecretFinder {
    // Simple exact fraction with BigInteger
    static class BigFraction {
        BigInteger num;
        BigInteger den;
        BigFraction(BigInteger n, BigInteger d) {
            if (d.signum() == 0) throw new ArithmeticException("zero denominator");
            if (d.signum() < 0) { n = n.negate(); d = d.negate(); }
            BigInteger g = n.gcd(d);
            if (!g.equals(BigInteger.ONE)) { n = n.divide(g); d = d.divide(g); }
            this.num = n; this.den = d;
        }
        BigFraction(BigInteger n) { this(n, BigInteger.ONE); }
        BigFraction add(BigFraction other) {
            BigInteger n = this.num.multiply(other.den).add(other.num.multiply(this.den));
            BigInteger d = this.den.multiply(other.den);
            return new BigFraction(n, d);
        }
        BigFraction multiply(BigFraction other) {
            return new BigFraction(this.num.multiply(other.num), this.den.multiply(other.den));
        }
        BigFraction divide(BigFraction other) {
            if (other.num.signum() == 0) throw new ArithmeticException("divide by zero");
            return new BigFraction(this.num.multiply(other.den), this.den.multiply(other.num));
        }
        public String toString() {
            if (den.equals(BigInteger.ONE)) return num.toString();
            return num.toString() + "/" + den.toString();
        }
    }

    // decode a string value in given base to BigInteger
    static BigInteger decodeValue(String s, int base) {
        s = s.trim().toLowerCase();
        BigInteger result = BigInteger.ZERO;
        BigInteger b = BigInteger.valueOf(base);
        for (char ch : s.toCharArray()) {
            int digit = Character.digit(ch, base);
            if (digit == -1) {
                throw new IllegalArgumentException("Invalid digit '" + ch + "' for base " + base);
            }
            result = result.multiply(b).add(BigInteger.valueOf(digit));
        }
        return result;
    }

    // compute f(0) using Lagrange interpolation (points list length = k)
    static BigFraction computeF0(List<Map.Entry<Integer, BigInteger>> pts) {
        int k = pts.size();
        BigFraction f0 = new BigFraction(BigInteger.ZERO);
        for (int j = 0; j < k; ++j) {
            BigInteger xj = BigInteger.valueOf(pts.get(j).getKey());
            BigInteger yj = pts.get(j).getValue();
            BigFraction lj = new BigFraction(BigInteger.ONE);
            for (int i = 0; i < k; ++i) {
                if (i == j) continue;
                BigInteger xi = BigInteger.valueOf(pts.get(i).getKey());
                BigInteger num = xi.negate();              // (0 - xi)
                BigInteger den = xj.subtract(xi);         // (xj - xi)
                lj = lj.multiply(new BigFraction(num, den));
            }
            BigFraction term = lj.multiply(new BigFraction(yj));
            f0 = f0.add(term);
        }
        return f0;
    }

    // given a JsonObject representing the test-case, decode points and compute the secret
    static BigInteger secretFromJson(JsonObject root) {
        JsonObject keys = root.getAsJsonObject("keys");
        if (keys == null) throw new IllegalArgumentException("missing keys object");
        int k = keys.get("k").getAsInt();

        // collect points (x -> decoded y)
        List<Map.Entry<Integer, BigInteger>> points = new ArrayList<>();
        for (Map.Entry<String, JsonElement> e : root.entrySet()) {
            String key = e.getKey();
            if (key.equals("keys")) continue;
            if (key == null || key.trim().isEmpty()) continue;
            JsonObject valObj = e.getValue().getAsJsonObject();
            String baseStr = valObj.get("base").getAsString();
            String valStr  = valObj.get("value").getAsString();
            int x = Integer.parseInt(key);
            int base = Integer.parseInt(baseStr);
            BigInteger y = decodeValue(valStr, base);
            points.add(new AbstractMap.SimpleEntry<>(x, y));
        }
        // sort by x
        points.sort(Comparator.comparingInt(Map.Entry::getKey));
        if (points.size() < k) throw new IllegalArgumentException("not enough points, have " + points.size() + " need " + k);

        // take first k points (any k points would do if consistent)
        List<Map.Entry<Integer, BigInteger>> chosen = points.subList(0, k);
        BigFraction f0 = computeF0(chosen);
        if (!f0.den.equals(BigInteger.ONE)) {
            // Ideally should be integer per constraints; return fraction numerator/denominator if not integer
            throw new ArithmeticException("Interpolated f(0) is not integer: " + f0.toString());
        }
        return f0.num;
    }

    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            System.out.println("Usage: java -cp .:gson.jar SecretFinder <test1.json> [<test2.json> ...]");
            return;
        }
        Gson gson = new Gson();
        for (String path : args) {
            String text = new String(Files.readAllBytes(Paths.get(path)));
            JsonObject root = gson.fromJson(text, JsonObject.class);
            try {
                BigInteger secret = secretFromJson(root);
                System.out.println("File: " + path + "  -> secret = " + secret.toString());
            } catch (Exception ex) {
                System.out.println("File: " + path + "  -> ERROR: " + ex.getMessage());
            }
        }
    }
}
