import javax.swing.JFrame;

/**
 * @author Christoph Riesinger (riesinge@in.tum.de)
 * @author Jürgen Bräckle (braeckle@in.tum.de)
 * @author Sebastian Rettenberger (rettenbs@in.tum.de)
 * @since November 02, 2012
 * @version 1.1
 * 
 *          This class just contains a main() method to use the FastMath class
 *          and to invoke the plotter.
 */
public class Test_FastInverse {

	/** Beispielwerte fuer IEEE Standard mit 32 Bits */
	private static int MAGIC_NUMBER = 0x5f3759df;

	private static int anzBitsExponent = 4;	//8
	private static int anzBitsMantisse = 8;	//24

	/**
	 * Uses the FastMath class and invokes the plotter. In a logarithmically
	 * scaled system, the exakt solutions of 1/sqrt(x) are shown in green, and
	 * the absolute errors of the Fast Inverse Square Root in red. Can be used
	 * to test and debug the own implementation of the fast inverse square root
	 * algorithm and to play while finding an optimal magic number.
	 * 
	 * @param args
	 *            args is ignored.
	 */
	public static void main(String[] args) {

		Gleitpunktzahl.setSizeExponent(anzBitsExponent);
		Gleitpunktzahl.setSizeMantisse(anzBitsMantisse);
		//FastMath.setMagic(MAGIC_NUMBER);

		int numOfSamplingPts = 1001;
		float[] xData = new float[numOfSamplingPts];
		float[] yData = new float[numOfSamplingPts];
		float x = 0.10f;

		/* calculate data to plot */
		for (int i = 0; i < numOfSamplingPts; i++) {
			xData[i] = x;
			Gleitpunktzahl y = new Gleitpunktzahl(x);
			yData[i] = (float) FastMath.absInvSqrtErr(y); //abs -> rel

			x *= Math.pow(100.0d, 1.0d / numOfSamplingPts);
		}

		/* initialize plotter */
		JFrame frame = new JFrame();
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		try {
			frame.add(new Plotter(xData, yData));
		} catch (InstantiationException exception) {
			exception.printStackTrace();
			System.exit(1);
		}
		frame.setSize(960, 720);
		frame.setLocation(0, 0);
		frame.setVisible(true);
	}
}
