package planeten;

import java.awt.Color;
import java.util.ArrayList;

import ode.Einschrittverfahren;
import ode.ODE;

/**
 * Diese Klasse verwaltet ein System von Planeten
 * 
 * @author braeckle
 * 
 */
public class Planetensystem {

	ArrayList<Planet> planeten = new ArrayList<Planet>();

	ODE ode;

	public Planetensystem() {
		addPlanet("Sonne", 1.9891E30, 1392700, new Color(224, 215, 75), 0, 0,
				0, 0);
		addPlanet("Merkur", 3.302E23, 4879.4, new Color(166, 115, 53),
				57.909E6, 0, 0, -47.87);
		addPlanet("Venus", 4.869E24, 12103.6, new Color(120, 74, 32), 108.16E6,
				0, 0, -35.02);
		addPlanet("Erde", 5.974E24, 12756.32, new Color(2, 60, 106), 149.6E6,
				0, 0, -29.78);
		addPlanet("Mars", 6.419E23, 6792.4, new Color(155, 129, 104), 227.99E6,
				0, 0, -24.13);
		addPlanet("Jupiter", 1.899E27, 142984, new Color(148, 144, 145),
				778.36E6, 0, 0, -13.07);
		addPlanet("Saturn", 5.685E26, 120536, new Color(179, 160, 129),
				1433.4E6, 0, 0, -9.69);
		addPlanet("Uranus", 8.683E25, 51118, new Color(192, 230, 233),
				2872.4E6, 0, 0, -6.81);
		addPlanet("Neptun", 1.0243E26, 49528, new Color(102, 155, 241), 4495E6,
				0, 0, -5.43);

		addPlanet("Mond", 7.349E22, 3476, new Color(1, 1, 0), 149.6E6, 384400,
				1.023, -29.78);
		addPlanet("Hale-Bopp", 1.3E16, 60, new Color(0, 0, 0),
				-0.914 * 149.6E6, 0, 0, 44.005);
		addPlanet("HalleyscherKomet", 2E14, 10, new Color(0, 0, 0),
				-0.568 * 149.6E6, 0, 0, 54.57);

		/* setzt die ODE */
		ode = new NewtonGravitationODE(planeten);
	}

	public void addPlanet(String name, double mass, double durchmesser,
			Color color, double initPosX, double initPosY, double initVelX,
			double initVelY) {
		Planet p = new Planet(name, mass, durchmesser, color);
		double[] initPos = new double[2];
		double[] initVel = new double[2];
		initPos[0] = initPosX;
		initPos[1] = initPosY;
		initVel[0] = initVelX;
		initVel[1] = initVelY;
		p.setPos(initPos);
		p.setVel(initVel);
		planeten.add(p);
	}

	/**
	 * fuehrt einen Zeitschritt des uebergebenen Verfahrens auf die
	 * NewtonGravitation-ODE aus
	 */
	public void nextStep(Einschrittverfahren method, double delta_t) {
		double[] y_k = new double[4 * planeten.size()];
		for (int i = 0; i < planeten.size(); i++) {
			Planet p = planeten.get(i);
			y_k[4 * i] = p.pos[0];
			y_k[4 * i + 1] = p.pos[1];
			y_k[4 * i + 2] = p.vel[0];
			y_k[4 * i + 3] = p.vel[1];
		}
		/* Einschrittverfahren anwenden */
		double[] y = method.nextStep(y_k, 0, delta_t, ode);

		for (int i = 0; i < planeten.size(); i++) {
			Planet p = planeten.get(i);
			p.pos[0] = y[4 * i];
			p.pos[1] = y[4 * i + 1];
			p.vel[0] = y[4 * i + 2];
			p.vel[1] = y[4 * i + 3];
		}
	}
}

/**
 * Diese ODE beschreibt die Bahnen der einzelnen Planeten, basierend auf den
 * Gravitationskraeften der Planeten untereinander. Die Länge von y ist ein
 * Vielfaches von 4. Jeweils 4 aufeinanderfolgende Werte beschreiben ein
 * Objekt(x-Koord, y-Koord, Geschwindigkeit in x- und y-Richtung)
 * 
 * Alle Zeiten t in Sekunden
 */
class NewtonGravitationODE implements ODE {

	/** Gravitationskonstante in km^3/kg s^2 */
	final static double G = 6.673E-20;

	ArrayList<Planet> planeten;

	public NewtonGravitationODE(ArrayList<Planet> planeten) {
		this.planeten = planeten;
	}

	@Override
	public double[] auswerten(double t, double[] y) {
		int n = y.length / 4;
		double[] yPunkt = new double[4 * n];
		for (int i = 0; i < n; i++) {
			yPunkt[4 * i] = y[4 * i + 2];
			yPunkt[4 * i + 1] = y[4 * i + 3];

			yPunkt[4 * i + 2] = 0;
			yPunkt[4 * i + 3] = 0;
			for (int j = 0; j < n; j++) {
				if (i != j) {
					double mj = planeten.get(j).mass;
					double distance = Math.sqrt(Math
							.pow(y[4 * i] - y[4 * j], 2)
							+ Math.pow(y[4 * i + 1] - y[4 * j + 1], 2));
					double factor = mj / Math.pow(distance, 3);
					yPunkt[4 * i + 2] += factor * (y[4 * i] - y[4 * j]);
					yPunkt[4 * i + 3] += factor * (y[4 * i + 1] - y[4 * j + 1]);
				}
			}
			yPunkt[4 * i + 2] *= -G;
			yPunkt[4 * i + 3] *= -G;
		}

		return yPunkt;
	}
}
