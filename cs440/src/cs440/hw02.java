package cs440;

public class hw02 {
	public static void main(String[] args) {

		toString(getNeuralNetsOutput(0, 1, 0));
		
		toString(getNeuralNetsOutput(1, 0, 0));
		

		
	}
	
	public static void toString(int[] a){
		for (int i=0; i<a.length-1; i++){
			System.out.print(a[i] + ", ");
		}
		System.out.println(a[a.length-1]);
	}
	
	
	
	public static int[] getNeuralNetsOutput(int a, int b, int c){
		int[] ans = new int[3];
		int[] Nodes = new int[5];
		if (a*1+b*1 > 0.5){
			Nodes[0] = 1;
		}
		if (b*1+c*1 >0.5){
			Nodes[1] =1;
		}
		Nodes[2] = Nodes[0]*1+Nodes[1]*-1;
		Nodes[3] = Nodes[0]*1+Nodes[1]*1;
		Nodes[4] = Nodes[0]*-1+Nodes[1]*1;
		if(Nodes[2]>0.5){
			Nodes[2] = 1;
		} else {
			Nodes[2]=0;
		}
		if(Nodes[3]>1.5){
			Nodes[3] = 1;
		} else {
			Nodes[3] = 0;
		}
		if(Nodes[4]>0.5){
			Nodes[4] = 1;
		} else {
			Nodes[4] = 0;
		}
		ans[0]=Nodes[2];
		ans[1]=Nodes[3];
		ans[2]=Nodes[4];
		return ans;
	}
}
