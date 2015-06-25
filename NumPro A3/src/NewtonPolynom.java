import java.util.Arrays;

/**
 * Die Klasse Newton-Polynom beschreibt die Newton-Interpolation. Die Klasse
 * bietet Methoden zur Erstellung und Auswertung eines Newton-Polynoms, welches
 * uebergebene Stuetzpunkte interpoliert.
 * 
 * @author braeckle
 * 
 */
public class NewtonPolynom implements InterpolationMethod {

	/** Stuetzstellen xi */
	double[] x;

	/**
	 * Koeffizienten/Gewichte des Newton Polynoms p(x) = a0 + a1*(x-x0) +
	 * a2*(x-x0)*(x-x1)+...
	 */
	double[] a;

	/**
	 * die Diagonalen des Dreiecksschemas. Diese dividierten Differenzen werden
	 * fuer die Erweiterung der Stuetzstellen benoetigt.
	 */
	double[] f;

	/**
	 * leerer Konstruktore
	 */
	public NewtonPolynom() {
	};

	/**
	 * Konstruktor
	 * 
	 * @param x
	 *            Stuetzstellen
	 * @param y
	 *            Stuetzwerte
	 */
	public NewtonPolynom(double[] x, double[] y) {
		this.init(x, y);
	}

	/**
	 * {@inheritDoc} Zusaetzlich werden die Koeffizienten fuer das
	 * Newton-Polynom berechnet.
	 */
	@Override
	public void init(double a, double b, int n, double[] y) {
		x = new double[n + 1];
		double h = (b - a) / n;

		for (int i = 0; i < n + 1; i++) {
			x[i] = a + i * h;
		}
		computeCoefficients(y);
	}

	/**
	 * Initialisierung der Newtoninterpolation mit beliebigen Stuetzstellen. Die
	 * Faelle "x und y sind unterschiedlich lang" oder "eines der beiden Arrays
	 * ist leer" werden nicht beachtet.
	 * 
	 * @param x
	 *            Stuetzstellen
	 * @param y
	 *            Stuetzwerte
	 */
	public void init(double[] x, double[] y) {
		this.x = Arrays.copyOf(x, x.length);
		computeCoefficients(y);
	}

	/**
	 * computeCoefficients belegt die Membervariablen a und f. Sie berechnet zu
	 * uebergebenen Stuetzwerten y, mit Hilfe des Dreiecksschemas der
	 * Newtoninterpolation, die Koeffizienten a_i des Newton-Polynoms. Die
	 * Berechnung des Dreiecksschemas soll dabei lokal in nur einem Array der
	 * Laenge n erfolgen (z.B. spaltenweise Berechnung). Am Ende steht die
	 * Diagonale des Dreiecksschemas in der Membervariable f, also f[0],f[1],
	 * ...,f[n] = [x0...x_n]f,[x1...x_n]f,...,[x_n]f. Diese koennen spaeter bei
	 * der Erweiterung der Stuetzstellen verwendet werden.
	 * 
	 * Es gilt immer: x und y sind gleich lang.
	 */
	private void computeCoefficients(double[] y) {
		
		//NOTE: as said in lecture, Newton is not suitable for anything but very small pictures
		
		//tested to be correct: passed
		int n = x.length;
		double[][] tri = new double[n][n];
		
		//for each line in the Dreiecksschema:
		//	step 1: copy yi into the first column
		for(int i = 0; i < n; i++) {
			tri[i][0] = y[i];
		}
		for(int k = 1; k<n; k++)
		{
			for(int i= 0; i<(n-k); i++)
			{
				tri[i][k] = (tri[i+1][k-1]-tri[i][k-1])/(x[i+k]-x[i]);
			}
		}
		//MatrixPlotter.plot(tri, n);
		f = new double[n];
		a = new double[n];
		//copy right/down diagonal into 'f' array, copy first line into 'a' array
		for(int i = 0; i < f.length; i++) {
			f[i] = tri[(f.length-1)-i][i];
			a[i] = tri[0][i];
		}
	}

	/**
	 * Gibt die Koeffizienten des Newton-Polynoms a zurueck
	 */
	public double[] getCoefficients() {
		return a;
	}

	/**
	 * Gibt die Dividierten Differenzen der Diagonalen des Dreiecksschemas f
	 * zurueck
	 */
	public double[] getDividedDifferences() {
		return f;
	}

	/**
	 * addSamplingPoint fuegt einen weiteren Stuetzpunkt (x_new, y_new) zu x
	 * hinzu. Daher werden die Membervariablen x, a und f vergoessert und
	 * aktualisiert . Das gesamte Dreiecksschema muss dazu nicht neu aufgebaut
	 * werden, da man den neuen Punkt unten anhaengen und das alte
	 * Dreiecksschema erweitern kann. Fuer diese Erweiterungen ist nur die
	 * Kenntnis der Stuetzstellen und der Diagonalen des Schemas, bzw. der
	 * Koeffizienten noetig. Ist x_new schon als Stuetzstelle vorhanden, werden
	 * die Stuetzstellen nicht erweitert.
	 * 
	 * @param x_new
	 *            neue Stuetzstelle
	 * @param y_new
	 *            neuer Stuetzwert
	 */
	public void addSamplingPoint(double x_new, double y_new) {
		
		/* Tested */
		int size = x.length;
		double[] newx = new double[size+1];
		
		for(int i = 0; i<size; i++) {
			if(x[i] == x_new)
				return;
			else
				newx[i]=x[i];
		}
		newx[size]= x_new;
		x = newx;
		
		double[] f_new = new double[f.length+1];
		double[] a_new = new double[a.length+1];
		f_new[0]=y_new;
		for(int i =1;i<f.length+1; i++)
		{
			f_new[i] = (f_new[i-1] - f[i-1]) / (x[0] - x[i]); 
			a_new[i] = a[i]; 
		}
		a_new[f.length+1] = f_new[f.length+1]; //loop goes till i = f.lenth
		f = f_new;
		a = a_new;
	}

	/**
	 * {@inheritDoc} Das Newton-Polynom soll effizient mit einer Vorgehensweise
	 * aehnlich dem Horner-Schema ausgewertet werden. Es wird davon ausgegangen,
	 * dass die Stuetzstellen nicht leer sind.
	 */
	@Override
	public double evaluate(double z) {
		
		/* Tested */
		double result=0;
		for(int i=0; i<f.length; i++)
		{
			double product = a[i];
			for(int j=0; j<i;j++)
			{
				product*=(z-x[j]);
			}
			result+=product;
		}
		return result;
	}
}
