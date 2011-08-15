import java.util.Random;


public class randtest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Random rand=new Random(0);
		for(int i = 0; i < 10; i++){
			System.out.println(rand.nextInt());
		}
	}

}
