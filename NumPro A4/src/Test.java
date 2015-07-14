import java.util.Arrays;

import ode.ExpliziterEuler;
import ode.Heun;
import ode.ODE;
import ode.RungeKutta4;
import planeten.PlanetenGUI;
import freierfall.FastTransportGui;
import dft.DFT;
import dft.IFFT;
import dft.Complex;

public class Test {

	/**
	 * Hier werden die GUIs fuer die Freie-Fall- und die
	 * Planetensystemsimulation gestartet, und einzelne Testfaelle
	 * durchgefuehrt.
	 */
	public static void main(String[] args) {

		/**************************************/
		boolean startPlanetensystem = true;
		boolean startFreierFall = true;

		boolean testExpliziteVerfahren = true;

		boolean testFFT = true;
		/**************************************/

		if (startPlanetensystem) {
			new PlanetenGUI().setVisible(true);
		}

		if (startFreierFall) {
			new FastTransportGui().setVisible(true);
		}

		if (testExpliziteVerfahren)
			testExpliziteVerfahren();

		if (testFFT)
			testFFT();
	}

	public static void testExpliziteVerfahren() {
		System.out.println("Teste Explizite Verfahren");
		ODE ode = new ODE() {

			@Override
			public double[] auswerten(double t, double[] y) {
				double[] v = new double[1];
				v[0] = t * y[0];
				return v;
			}
		};

		double delta_t = 1;
		double t0 = 0;
		double[] y0 = { 1 };

		/* Expl Euler */
		System.out.println("\nTeste Expliziten Euler.");
		ExpliziterEuler euler = new ExpliziterEuler();
		double[] y = Arrays.copyOf(y0, y0.length);
		double t = t0;
		for (int k = 1; k <= 4; k++) {
			y = euler.nextStep(y, t, delta_t, ode);
			System.out.println("y" + k + " = " + y[0]);
			t = t + delta_t;
		}
		System.out.println("Richtig waere: Eigene Beispiele ueberlegen");

		/* Heun */
		System.out.println("\nTeste Heun.");
		Heun heun = new Heun();
		y = Arrays.copyOf(y0, y0.length);
		t = t0;
		for (int k = 1; k <= 4; k++) {
			y = heun.nextStep(y, t, delta_t, ode);
			System.out.println("y" + k + " = " + y[0]);
			t = t + delta_t;
		}
		System.out.println("Richtig waere: Eigene Beispiele ueberlegen");

		/* Runge Kutta4 */
		System.out.println("\nTeste Runge-Kutta4.");
		RungeKutta4 rk4 = new RungeKutta4();
		y = Arrays.copyOf(y0, y0.length);
		t = t0;
		for (int k = 1; k <= 4; k++) {
			y = rk4.nextStep(y, t, delta_t, ode);
			System.out.println("y" + k + " = " + y[0]);
			t = t + delta_t;
		}
		System.out
				.println("Richtig waeren gerundet: Eigene Beispiele ueberlegen");

		System.out.println("*************************************\n");
	}

	public static void testFFT() {
		System.out.println("Teste Fast Fourier Transformation");

		double[] v = new double[4];
		for (int i = 0; i < 4; i++)
			v[i] = i + 1;
		Complex[] c = dft.DFT.dft(v);
		Complex[] v2 = dft.IFFT.ifft(c);

		for (int i = 0; i < 4; i++) {
			System.out.println(v2[i]);
			System.out.println(v[i]);
		}
		System.out
				.println("Richtig waeren gerundet: Eigene Beispiele ueberlegen");

		System.out.println("*************************************\n");
	}

}
