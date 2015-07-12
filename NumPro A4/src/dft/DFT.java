package dft;

/**
 * Diskrete Fourier-Transformation
 *
 * @author Sebastian Rettenberger
 */
public class DFT {
	/**
	 * Diskrete Fourier-Transformation (DFT)
	 */
	public static Complex[] dft(double[] v) {
		int n = v.length;
		Complex[] c = new Complex[n];

		Complex omega = Complex.fromPolar(1, 2 * Math.PI / n).conjugate();

		for (int i = 0; i < n; i++) {
			c[i] = new Complex();

			for (int j = 0; j < n; j++) {
				Complex t = omega.power(i * j);
				t = t.mul(new Complex(v[j] / n));

				c[i] = c[i].add(t);
			}
		}

		return c;
	}
}