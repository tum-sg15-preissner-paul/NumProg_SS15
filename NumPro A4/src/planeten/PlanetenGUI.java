package planeten;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Ellipse2D;
import java.util.Timer;
import java.util.TimerTask;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuBar;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

import ode.Einschrittverfahren;
import ode.ExpliziterEuler;
import ode.Heun;
import ode.RungeKutta4;

/**
 * Diese Benutzeroberflaeche zeigt die Simulation eines Planetensystems
 * 
 * @author braeckle
 * 
 */

/**
 * ViewComponent verwaltet die Anzeigeflaeche fuer das Planetensystem
 */
class ViewComponent extends JComponent {
	private static final long serialVersionUID = 1L;

	/** Die angezeigten Planeten */
	public Planetensystem planetensystem;

	/* offset versetzt das Feld (0,0) in die Mitte des Bildes */
	int[] offset = new int[2];

	/* groesster Durchmesser im System */
	double largestDiameter;
	/* groesster Durchmesser im System soll dieser Pixelzahl entsprechen */
	int largestDiameterInPixel = 50;
	/* Abstand zweier Pixel entspricht pixellaenge km */
	double pixellaenge; // km/pixel

	/* Zeit, die seit Start der Simulation vergangen ist, in sec */
	double timecounter = 0;

	public ViewComponent() {
		super();
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		offset[0] = Math.min(dim.width * 9 / 20, 400);
		offset[1] = Math.min(dim.height * 9 / 20, 320);
	}

	public void setPlanetensystem(Planetensystem p) {
		timecounter = 0;
		planetensystem = p;

		/* finde groesste Entfernung vom Zentrum */
		double maxDist = 0;
		for (Planet planet : planetensystem.planeten) {
			double dist = Math.sqrt(Math.pow(planet.pos[0], 2)
					+ Math.pow(planet.pos[1], 2));
			if (dist > maxDist)
				maxDist = dist;
		}
		/*
		 * Setzte entsprechend dieser Entfernung die die Maße des gezeigten
		 * Ausschnitts
		 */
		pixellaenge = Math.min(offset[0], offset[1]) / maxDist;

		/* finde groessten Planeten */
		double maxdur = 0;
		for (Planet planet : planetensystem.planeten) {
			double dur = planet.durchmesser;
			if (dur > maxdur)
				maxdur = dur;
		}
		largestDiameter = maxdur;
	}

	public void nextStep(Einschrittverfahren method, double timestep) {
		timecounter += timestep;
		planetensystem.nextStep(method, timestep);
		repaint();
	}

	@Override
	protected void paintComponent(Graphics g) {
		Graphics2D g2d = (Graphics2D) g;
		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
				RenderingHints.VALUE_ANTIALIAS_ON);
		/* Zeitangabe */
		double days = timecounter / 86400;
		int years = (int) days / 365;
		char[] data = (years + " Jahre, " + (int) days % 365 + " Tage")
				.toCharArray();
		g2d.drawChars(data, 0, data.length, 100, 100);
		/* Planeten zeichnen */
		for (Planet planet : planetensystem.planeten) {
			g2d.setColor(planet.color);
			// double durchmesser = Math.log(planet.durchmesser) /
			// Math.log(largestInKM) * largestInPixel;
			// double durchmesser = Math.sqrt(planet.durchmesser) /
			// Math.sqrt(largestInKM) * largestInPixel;
			double durchmesser = Math.sqrt(planet.durchmesser)
					/ Math.sqrt(largestDiameter)
					/ (Math.log(planet.durchmesser) / Math.log(largestDiameter))
					* largestDiameterInPixel + 2;
			double x = planet.pos[0] * pixellaenge + offset[0] - durchmesser
					/ 2;
			double y = planet.pos[1] * pixellaenge + offset[1] - durchmesser
					/ 2;
			Ellipse2D.Double circle = new Ellipse2D.Double(x, y, durchmesser,
					durchmesser);
			g2d.fill(circle);
		}
	}
}

/***********************************************************/

/** Die Benutzeroberflaeche */
public class PlanetenGUI extends JFrame implements ActionListener {
	private static final long serialVersionUID = 1L;

	private ViewComponent viewComponent = new ViewComponent();

	Timer timer = new Timer();
	boolean simulationlaeuft = false;

	JMenuBar mbar = new JMenuBar();

	/** Combobox fuer die Einschrittmethoden */
	String[] einschrittStrings = { "ExpliziterEuler", "Heun", "RungaKutta" };
	JComboBox einschrittBox = new JComboBox(einschrittStrings);

	/** Textfeld fuer die Eingabe des Zeitschritts in Tagen */
	JTextField timestepField = new JTextField(1 + "", 4);

	/** Textfeld fuer die Eingabe der Pause zwischen zwei Zeitschritten */
	JTextField delayField = new JTextField(10 + "", 4);

	/** Der Start-Stop-Button */
	JButton startStopButton = new JButton("Start");

	/** Der Reset-Button */
	JButton resetButton = new JButton("Reset");

	/**
	 * Der Konstruktor fuegt alle Komponenten zusammen
	 */
	public PlanetenGUI() {
		super("Planetensystem");

		viewComponent.setPlanetensystem(new Planetensystem());

		GridLayout layout = new GridLayout(3, 3);

		mbar.setLayout(layout);

		einschrittBox.setSelectedIndex(0);
		startStopButton.addActionListener(this);
		resetButton.addActionListener(this);

		mbar.add(startStopButton);
		mbar.add(einschrittBox);
		mbar.add(resetButton);

		mbar.add(new JLabel("Zeitschritt: ", SwingConstants.RIGHT));
		mbar.add(timestepField);
		mbar.add(new JLabel("Tage", SwingConstants.LEFT));

		mbar.add(new JLabel("Delay: ", SwingConstants.RIGHT));
		mbar.add(delayField, BorderLayout.SOUTH);
		mbar.add(new JLabel("msecs", SwingConstants.LEFT));

		setJMenuBar(mbar);
		add(viewComponent);

		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		setSize((int) Math.min(dim.width * 9 / 10, 800),
				(int) Math.min(dim.height * 9 / 10, 640));
	}

	public void actionPerformed(ActionEvent e) {

		JComponent source = (JComponent) (e.getSource());

		if (source == (JComponent) resetButton) {
			timer.cancel();
			viewComponent.setPlanetensystem(new Planetensystem());
			viewComponent.repaint();
			simulationlaeuft = false;
			startStopButton.setText("Start");

		} else if (source == (JComponent) startStopButton) {

			if (!simulationlaeuft) {
				timer = new Timer();
				simulationlaeuft = true;
				startStopButton.setText("Stop");

				TimerTask task = new TimerTask() {

					@Override
					public void run() {
						if (!simulationlaeuft)
							timer.cancel();
						else {
							int mode = einschrittBox.getSelectedIndex();
							Einschrittverfahren method;
							switch (mode) {
							case 0:
								method = new ExpliziterEuler();
								break;
							case 1:
								method = new Heun();
								break;
							case 2:
								method = new RungeKutta4();
								break;
							default:
								method = new ExpliziterEuler();
							}
							double timestep;
							try {
								timestep = Double.parseDouble(timestepField
										.getText()) * 86400;
							} catch (NumberFormatException e) {
								timestep = 86400;
							}
							viewComponent.nextStep(method, timestep);
						}
					}
				};

				int delay;
				try {
					delay = Integer.parseInt(delayField.getText());
					if (delay < 1) {
						delay = 10;
						this.delayField.setText("" + delay);
					}
				} catch (NumberFormatException f) {
					delay = 10;
					this.delayField.setText("" + delay);
				}

				timer.schedule(task, 0, (long) delay);
			} else {
				// SToppe Simulation
				simulationlaeuft = false;
				startStopButton.setText("Start");
			}
		}

	}

	/*********************************************************/
	/** Starten der Benutzeroberflaeche */
	public static void main(String[] args) {
		new PlanetenGUI().setVisible(true);
	}

}