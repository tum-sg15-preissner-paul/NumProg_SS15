package ode;

import java.util.Arrays;

/**
 * Das Einschrittverfahren "Expliziter Euler"
 * 
 * @author braeckle
 * 
 */
public class ExpliziterEuler implements Einschrittverfahren {

	public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) {
		//tests suggest it works, + Tutor exercise provides same solution
		
		//f = f(t,r(t)) as in task sheet
		double[] f = ode.auswerten(t, y_k);
		for(int i = 0; i < y_k.length; i++) {
			y_k[i] = y_k[i] + (delta_t * f[i]);			
		}
		
		return Arrays.copyOf(y_k, y_k.length);
	}

}
