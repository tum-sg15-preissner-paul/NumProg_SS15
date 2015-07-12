package ode;

import java.util.Arrays;

/**
 * Der klassische Runge-Kutta der Ordnung 4
 * 
 * @author braeckle
 * 
 */
public class RungeKutta4 implements Einschrittverfahren {

	@Override
	/**
	 * {@inheritDoc}
	 * Bei der Umsetzung koennen die Methoden addVectors und multScalar benutzt werden.
	 */
	public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) {
		// TODO: test
		
		//vectors needed for intermediate steps
		double[] k1 = new double[y_k.length];
		double[] k2 = new double[y_k.length];
		double[] k3 = new double[y_k.length];
		double[] k4 = new double[y_k.length];
		double[] f_k1, f_k2, f_k3, f_k4;
		
		//calc k1
		f_k1 = ode.auswerten(t, y_k);
		for(int i = 0; i < y_k.length; i++) {
			k1[i] = delta_t * f_k1[i];
		}
		
		//calc k2
		f_k2 = ode.auswerten(t + delta_t/2, addVectors(y_k, multScalar(k1, 0.5d)));
		for(int i = 0; i < y_k.length; i++) {
			k2[i] = delta_t * f_k2[i];
		}
		
		//calc k3
		f_k3 = ode.auswerten(t + delta_t/2, addVectors(y_k, multScalar(k2, 0.5d)));
		for(int i = 0; i < y_k.length; i++) {
			k3[i] = delta_t * f_k3[i];
		}
		
		//calc k4
		f_k4 = ode.auswerten(t + delta_t, addVectors(y_k, k3));
		for(int i = 0; i < y_k.length; i++) {
			k4[i] = delta_t * f_k4[i];
		}
		
		//calc new position
		for(int i = 0; i < y_k.length; i++) {
			y_k[i] = y_k[i] + (1d/6d) *
					addVectors(k1, //k1 +
							addVectors(multScalar(k2, 2d), //2*k2 +
									addVectors(multScalar(k3, 2d), k4)))[i]; //2*k3 + k4
		}
		
		return Arrays.copyOf(y_k, y_k.length);
	}

	/**
	 * addiert die zwei Vektoren a und b
	 */
	private double[] addVectors(double[] a, double[] b) {
		double[] erg = new double[a.length];
		for (int i = 0; i < a.length; i++) {
			erg[i] = a[i] + b[i];
		}
		return erg;
	}

	/**
	 * multipliziert den Skalar scalar auf den Vektor a
	 */
	private double[] multScalar(double[] a, double scalar) {
		double[] erg = new double[a.length];
		for (int i = 0; i < a.length; i++) {
			erg[i] = scalar * a[i];
		}
		return erg;
	}

}
