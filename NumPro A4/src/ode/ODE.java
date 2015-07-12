package ode;

/**
 * Diese Klasse beschreibt eine gewoehnliche Differentialgleichung y' = f(t,y),
 * deren rechte Seite ausgewertet werden kann.
 * 
 * @author braeckle
 * 
 */
public interface ODE {
	/**
	 * wertet die rechte Seite der ODE zu einem Zeitpunkt t mit einer aktuellen
	 * Loesung y aus
	 */
	public double[] auswerten(double t, double[] y);
}
