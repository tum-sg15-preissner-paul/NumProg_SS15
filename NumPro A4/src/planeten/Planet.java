package planeten;

import java.awt.Color;
import java.util.Arrays;

/** diese Klasse beschreibt einen Planeten */

public class Planet {

	String name;
	double mass; // in kg
	double durchmesser; // in km
	Color color;

	/** Aktuelle Position und Geschwindigket in x- und y-Richtung */
	double[] pos = new double[2]; // in km
	double[] vel = new double[2]; // in km/s

	public Planet(String name, double mass, double durchmesser, Color color) {
		this.name = name;
		this.mass = mass;
		this.color = color;
		this.durchmesser = durchmesser;
	}

	public void setPos(double[] pos) {
		if (pos.length == 2)
			this.pos = Arrays.copyOf(pos, 2);
	}

	public void setVel(double[] vel) {
		if (vel.length == 2)
			this.vel = Arrays.copyOf(vel, 2);
	}
}
