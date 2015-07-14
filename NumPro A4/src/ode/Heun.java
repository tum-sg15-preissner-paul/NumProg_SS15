package ode;

import java.util.Arrays;

/**
 * Das Einschrittverfahren von Heun
 * 
 * @author braeckle
 * 
 */
public class Heun implements Einschrittverfahren {

	@Override
	/**
	 * {@inheritDoc} 
	 * Nutzen Sie dabei geschickt den Expliziten Euler.
	 */
	public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) {
		//tests suggest it works, + Tutor exercise provides same solution
		
		//f = f(t,r(t)) as in task sheet
		double[] f = ode.auswerten(t, y_k);
		//r_ = r(t+deltaT)_ = r(t) + deltaT*f(t,r(t))
		double[] r_ = new double[y_k.length];
		
		//calculate first Heun step, r(t+deltaT)_
		for(int i = 0; i < y_k.length; i++) {
			r_[i] = y_k[i] + delta_t * f[i];
		}
		
		//f_delt = f(t+deltaT, r_)
		double[] f_delt = ode.auswerten(t + delta_t, r_);
		for(int i = 0; i < y_k.length; i++) {
			y_k[i] = y_k[i] + (delta_t/2) * (f[i] + f_delt[i]);
		}
		
		return Arrays.copyOf(y_k, y_k.length);
	}

}
