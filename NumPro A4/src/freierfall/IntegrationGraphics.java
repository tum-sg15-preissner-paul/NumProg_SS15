package freierfall;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Point;

import javax.swing.BorderFactory;
import javax.swing.JPanel;

import ode.ExpliziterEuler;
import ode.Heun;
import ode.ODE;
import ode.RungeKutta4;

/**
 * this class is responsible for drawing the ODE functions
 */
class IntegrationGraphics extends JPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	// / timestep size
	double timestepSize = 0.1;

	// / start position to start simulation
	double[] startPosition = new double[2];

	// / start velocity to start simulation
	double[] startVelocity = new double[2];

	// integration method to compute values for and to visualize
	boolean eulerIntegration = true;
	boolean heunIntegration = true;
	boolean rungeKuttaIntegration = true;

	// / computed integration points which are going to be drawn on screen
	double[][][] integrationPointsAll = new double[5][10000][2];
	int[] numberOfValidIntegrationPointsAll = new int[5];

	// / recently computed point of hit (linear interpolated position when
	// crossing the X-axis)
	double[] hitPointsX = new double[5];

	// / some drawing properties
	Dimension viewportSize = new Dimension(800, 600);
	Point origin = new Point(viewportSize.width / 2, 50);

	// zu loesende ODE
	final double gravity = 9.81;
	int odeNr;
	ODE ode;

	// / Graphics handler
	Graphics graphics;

	/**
	 * constructor
	 */
	public IntegrationGraphics() {
		super();
		setOpaque(false); // we don't paint all our bits
		setLayout(new BorderLayout());
		setBorder(BorderFactory.createLineBorder(Color.black));
	}

	/**
	 * we prefer a drawing area of 800x600
	 */
	public Dimension getPreferredSize() {
		return new Dimension(800, 600);
	}

	/**
	 * update the positions for given integration method 0: runge kutta with
	 * fine timesteps 1: euler 2: heun 3: runge-kutta
	 */
	public void computeIntegratorValues(int integrationMethod) {

		updateODE();

		double[] pos = startPosition.clone();

		double[][] integrationPoints = integrationPointsAll[integrationMethod];

		integrationPoints[0] = pos.clone();

		hitPointsX[integrationMethod] = -666.0;

		int i;
		boolean crossingFound = false;
		for (i = 1; i < 10000; i++) {
			double prevTime = (double) (i - 1) * timestepSize;

			switch (integrationMethod) {
			case 0:
				pos = new RungeKutta4().nextStep(pos, (double) (i - 1) * 0.01,
						0.01, ode);
				break;

			case 1:
				pos = new ExpliziterEuler().nextStep(pos, prevTime,
						timestepSize, ode);
				break;

			case 2:
				pos = new Heun().nextStep(pos, prevTime, timestepSize, ode);
				break;

			case 3:
				pos = new RungeKutta4().nextStep(pos, prevTime, timestepSize,
						ode);
				break;

			default:
				assert (false);
			}

			integrationPoints[i] = pos.clone();

			if (!crossingFound) {
				if (integrationPoints[i][1] <= 0.0) {
					crossingFound = true;

					double[] p1 = integrationPoints[i - 1];
					double[] p2 = integrationPoints[i];

					double deltay = (p1[1] - p2[1]);
					if (Math.abs(deltay) < 0.00001) {
						hitPointsX[integrationMethod] = p1[0];
					} else {
						hitPointsX[integrationMethod] = p1[0] + (p2[0] - p1[0])
								* (p1[1] / deltay);
					}
				}
			}
		}

		numberOfValidIntegrationPointsAll[integrationMethod] = i;
	}

	public void updateODE() {

		switch (odeNr) {
		case 0:
			ode = new ODE() {
				@Override
				public double[] auswerten(double t, double[] pos) {
					double[] v = new double[2];
					v[0] = startVelocity[0];
					v[1] = startVelocity[1] - gravity * t;
					return v;
				}
			};
			break;
		case 1:
			ode = new ODE() {
				@Override
				public double[] auswerten(double t, double[] pos) {
					double[] v = new double[2];
					v[0] = startVelocity[0] + pos[0] * pos[0] * 0.005;
					v[1] = startVelocity[1] - gravity * t - pos[0] * pos[0] * t
							* t * 0.001;
					return v;
				}
			};
			break;
		case 2:
			ode = new ODE() {
				@Override
				public double[] auswerten(double time, double[] pos) {
					double[] v = new double[2];
					v[0] = startVelocity[0] + time * 4.0;
					v[1] = startVelocity[1] + time * time * time * time * 0.1
							- time * time * 8.0;
					return v;
				}
			};
			break;
		case 3:
			ode = new ODE() {
				@Override
				public double[] auswerten(double t, double[] y) {
					double[] v = new double[2];
					v[0] = y[1] - 1;
					v[1] = -100 * (y[0] - 200) - 101 * y[1] + 1;
					return v;
				}
			};
			break;
		default:
			break;
		}
	}

	/**
	 * draw computed ODE values
	 */
	public void drawIntegratorValues(Color c, // /< color to draw function
			int integrationMethod, // /< integration method
			boolean drawCross // /< draw crosses to visualize discretized
								// timesteps
	) {
		double[][] integrationPoints = integrationPointsAll[integrationMethod];
		int numberOfValidIntegrationPoints = numberOfValidIntegrationPointsAll[integrationMethod];

		for (int i = 1; i < numberOfValidIntegrationPoints; i++) {
			double[] p1 = integrationPoints[i - 1];
			double[] p2 = integrationPoints[i];

			if (Math.abs(p2[0]) > 10000000)
				break;
			if (Math.abs(p2[1]) > 10000000)
				break;

			graphics.setColor(c);

			graphics.drawLine(origin.x + (int) p1[0], viewportSize.height
					- (origin.y + (int) p1[1]), origin.x + (int) p2[0],
					viewportSize.height - (origin.y + (int) p2[1]));

			if (drawCross) {
				graphics.setColor(Color.black);
				drawCross((int) p2[0], (int) p2[1], 2);
			}
		}
	}

	/**
	 * this method is executed to update the integrator values
	 */
	public void updateIntegratorValues() {
		computeIntegratorValues(0);

		if (eulerIntegration)
			computeIntegratorValues(1);

		if (heunIntegration)
			computeIntegratorValues(2);

		if (rungeKuttaIntegration)
			computeIntegratorValues(3);
	}

	/**
	 * callback handler from JPanel to update the graphics
	 */
	public void paintComponent(Graphics g) {
		graphics = g;

		// fill background
		g.setColor(Color.white);
		g.fillRect(0, 0, viewportSize.width, viewportSize.height);

		// draw axes
		g.setColor(Color.darkGray);
		g.drawLine(5, viewportSize.height - origin.y, viewportSize.width - 5,
				viewportSize.height - origin.y);
		g.drawLine(origin.x, 5, origin.x, viewportSize.height - 5);

		g.setColor(Color.red);
		g.fillOval(origin.x + (int) startPosition[0] - 5, viewportSize.height
				- (origin.y + (int) startPosition[1] + 5), 10, 10);

		drawIntegratorValues(Color.gray, 0, false);

		if (eulerIntegration)
			drawIntegratorValues(Color.red, 1, true);

		if (heunIntegration)
			drawIntegratorValues(Color.blue, 2, true);

		if (rungeKuttaIntegration)
			drawIntegratorValues(Color.green, 3, true);
	}

	/**
	 * draw a fancy cross with the given parameters
	 */
	public void drawCross(int x, // /< x-coordinate of center of cross
			int y, // /< y-coordinate of center of cross
			int size // /< size of cross
	) {
		graphics.drawLine(origin.x + x - size, viewportSize.height
				- (origin.y + y), origin.x + x + size, viewportSize.height
				- (origin.y + y));
		graphics.drawLine(origin.x + x, viewportSize.height
				- (origin.y + y + size), origin.x + x, viewportSize.height
				- (origin.y + y - size));
	}
}
