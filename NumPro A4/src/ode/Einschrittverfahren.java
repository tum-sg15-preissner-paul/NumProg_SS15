package ode;

/**
 * Das Interface fuer alle Einschrittverfahren.
 * 
 * @author braeckle
 * 
 */
public interface Einschrittverfahren {
	/**
	 * Diese Methode berechnet einen Schritt des Einschrittverfahrens
	 * 
	 * @param y_k
	 *            letzter berechneter Wert zum Zeitpunkt t (€ R^n)
	 * @param t
	 *            Zeitpunkt des letzten berechneten Schrittes
	 * @param delta_t
	 *            Zeitschritt zum naechsten Zeitpunkt
	 * @param ode
	 *            zu loesende ODE
	 * @return y_(k+1), Auswertung zum Zeitpunkt (t+delta_t)
	 */
	public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode);

}
