import java.util.Arrays;

/**
 * Die Klasse CubicSpline bietet eine Implementierung der kubischen Splines. Sie
 * dient uns zur effizienten Interpolation von aequidistanten Stuetzpunkten.
 * 
 * @author braeckle
 * 
 */
public class CubicSpline implements InterpolationMethod {

	/** linke und rechte Intervallgrenze x[0] bzw. x[n] */
	double a, b;

	/** Anzahl an Intervallen */
	int n;

	/** Intervallbreite */
	double h;

	/** Stuetzwerte an den aequidistanten Stuetzstellen */
	double[] y;

	/** zu berechnende Ableitungen an den Stuetzstellen */
	double yprime[];

	/**
	 * {@inheritDoc} Zusaetzlich werden die Ableitungen der stueckweisen
	 * Polynome an den Stuetzstellen berechnet. Als Randbedingungen setzten wir
	 * die Ableitungen an den Stellen x[0] und x[n] = 0.
	 */
	@Override
	public void init(double a, double b, int n, double[] y) {
		this.a = a;
		this.b = b;
		this.n = n;
		h = ((double) b - a) / (n);

		this.y = Arrays.copyOf(y, n + 1);

		/* Randbedingungen setzten */
		yprime = new double[n + 1];
		yprime[0] = 0;
		yprime[n] = 0;

		/* Ableitungen berechnen. Nur noetig, wenn n > 1 */
		if (n > 1) {
			computeDerivatives();
		}
	}

	/**
	 * getDerivatives gibt die Ableitungen yprime zurueck
	 */
	public double[] getDerivatives() {
		return yprime;
	}

	/**
	 * Setzt die Ableitungen an den Raendern x[0] und x[n] neu auf yprime0 bzw.
	 * yprimen. Anschliessend werden alle Ableitungen aktualisiert.
	 */
	public void setBoundaryConditions(double yprime0, double yprimen) {
		yprime[0] = yprime0;
		yprime[n] = yprimen;
		if (n > 1) {
			computeDerivatives();
		}
	}

	/**
	 * Berechnet die Ableitungen der stueckweisen kubischen Polynome an den
	 * einzelnen Stuetzstellen. Dazu wird ein lineares System Ax=c mit einer
	 * Tridiagonalen Matrix A und der rechten Seite c aufgebaut und geloest.
	 * Anschliessend sind die berechneten Ableitungen y1' bis yn-1' in der
	 * Membervariable yprime gespeichert.
	 * 
	 * Zum Zeitpunkt des Aufrufs stehen die Randbedingungen in yprime[0] und yprime[n].
	 * Speziell bei den "kleinen" Faellen mit Intervallzahlen n = 2
	 * oder 3 muss auf die Struktur des Gleichungssystems geachtet werden. Der
	 * Fall n = 1 wird hier nicht beachtet, da dann keine weiteren Ableitungen
	 * berechnet werden muessen.
	 */
	public void computeDerivatives() {
		/* TODO: test. fix if necessary */
		
		int nM = n-2;
		double[][] A = new double[nM][nM];
		double[] c = new double[nM];
		
		for(int i = 0; i < nM; i++) {
			/*-----Matrix A part*/
			A[i][i] = 4.0; //write 4 in diagonal A[i][i]
			if(i < nM-1)
				A[i+1][i] = A[i][i+1] = 1.0; //write 1 one to the right  and 1 one down from the 4
			/*-----*/
			
			/*-----Vector c part*/
			c[i]= y[i+2] - y[i];
		}
		c[0] -= h/3.0 * yprime[0];
		c[nM-1] -= h/3.0 * yprime[n];
			/*-----*/
		
		//Thomas-Algorithm, almost as seen on https://de.wikipedia.org/wiki/Thomas-Algorithmus
		/*NOTE: Only works when there's at least two fields to calculate?*/
		//TODO: Something is still off here, no idea what.
		if(n > 3) {
			//forward run
			double[] c_ = new double[nM-1];
			double[] d_ = new double[nM];
			for(int i = 0; i < nM; i++) {
				/*-----c'_i coefficients part*/
				if(i == 0) {
					c_[i] = 
							A[i][i+1] 
									/ A[i][i];
				} else if(i < nM-1) {
					c_[i] = A[i][i+1] / (A[i][i] - c_[i-1]*A[i][i-1]);
				}
				/*-----*/
				
				/*-----d'_i coefficients part*/
				if(i == 0) {
					d_[i] = c[i] / A[i][i];
				} else {
					c_[i] = (c[i] - c_[i-1]*A[i][i-1]) / (A[i][i] - c_[i-1]*A[i][i-1]);
				}
				/*-----*/
			}
			
			//backward run
			/*	x_n = d'_n,
				x_i = d'_i - (c'_i * x_{i+1}); i = n-1, n-2, ..., 1*/
			for(int i = nM; i > 0; i--) {
				//yprime's empty fields are in its range 1 to n-1, the i here is from 0 to n-2, so need to map accordingly
				if(i == nM) {
					yprime[i+1] = d_[i];
				} else {
					yprime[i+1] = d_[i] - (c_[i] * yprime[i+2]);
				}
			}
		} else if(n == 3){
			//n == 3, then only one or no yprime is left to solve?
			yprime[1] = (3.0/h * (y[2] - y[0] - (h/3.0)*yprime[0])) / 4.0;
		}
	}

	/**
	 * {@inheritDoc} Liegt z ausserhalb der Stuetzgrenzen, werden die
	 * aeussersten Werte y[0] bzw. y[n] zurueckgegeben. Liegt z zwischen den
	 * Stuetzstellen x_i und x_i+1, wird z in das Intervall [0,1] transformiert
	 * und das entsprechende kubische Hermite-Polynom ausgewertet.
	 */
	@Override
	public double evaluate(double z) {
		/* TODO: test. fix if necessary */
		
		if(z < a)
			return y[0];
		else if(z > b)
			return y[n];
		
		double i = Math.floor(z/h);
		z = (z-a)/(b-a); //assumption: Interval entirely positive
		
		//TODO: evaluate according Hermite Polynome
		//foo
		
		return 0.0;
	}
}
