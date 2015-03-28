package p03;

import java.util.Arrays;
import java.util.OptionalDouble;
import java.util.OptionalInt;

public class Testing {
	public static void main(String[] args) {
        int[] numbers = {1,2,3};
        OptionalInt highest = Arrays.stream(numbers).max();
        System.out.println("highest: " + highest);
        
        
        double[] a_ij_times_simga = new double[3];
        OptionalDouble highest1 = Arrays.stream(a_ij_times_simga).max();
        ;
	}

}
