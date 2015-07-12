package freierfall;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * GUI for visualization
 */
public class FastTransportPanel extends JPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	/**
	 * timestep size
	 */
	JSpinner spinnerTimestepSize;
	SpinnerNumberModel modelSpinnerTimestepSize;

	/**
	 * initial velocity in x direction
	 */
	JSpinner spinnerVelocityX;
	SpinnerNumberModel modelSpinnerVelocityX;

	/**
	 * initial velocity in y direction
	 */
	JSpinner spinnerVelocityY;
	SpinnerNumberModel modelSpinnerVelocityY;

	/**
	 * initial x coordinate of position
	 */
	JSpinner spinnerStartPositionX;
	SpinnerNumberModel modelSpinnerStartPositionX;

	/**
	 * initial y coordinate of position
	 */
	JSpinner spinnerStartPositionY;
	SpinnerNumberModel modelSpinnerStartPositionY;

	/**
	 * visualize / compute euler steps
	 */
	JCheckBox checkBoxEuler;
	JTextArea textAreaEuler;

	/**
	 * visualize / compute heun steps
	 */
	JCheckBox checkBoxHeun;
	JTextArea textAreaHeun;

	/**
	 * visualize / compute runge kutta steps
	 */
	JCheckBox checkBoxRungeKutta;
	JTextArea textAreaRungeKutta;

	IntegrationGraphics integrationGraphics = new IntegrationGraphics();

	/**
	 * constructor
	 */
	public FastTransportPanel() {
		super();

		// / initialize the panel gui components
		initPanel();
	}

	/**
	 * update values and repaint whole scene
	 */
	void updateGraphics() {
		integrationGraphics.updateIntegratorValues();

		textAreaEuler.setText("Euler X-zero-crossing: "
				+ String.valueOf((float) integrationGraphics.hitPointsX[1]));
		textAreaHeun.setText("Heun X-zero-crossing: "
				+ String.valueOf((float) integrationGraphics.hitPointsX[2]));
		textAreaRungeKutta.setText("Runge-Kutta4 X-zero-crossing: "
				+ String.valueOf((float) integrationGraphics.hitPointsX[3]));

		integrationGraphics.repaint();
	}

	/**
	 * initialize the panel with all GUI components
	 */
	public void initPanel() {
		JPanel rightPanel = new JPanel();

		GridBagConstraints c = new GridBagConstraints();
		c.gridy = 0;
		c.fill = GridBagConstraints.HORIZONTAL;

		rightPanel.setLayout(new GridBagLayout());

		/*
		 * FUNCTIONS
		 */
		ButtonGroup group = new ButtonGroup();

		JRadioButton fun0RadioButton = new JRadioButton("f(t,x) = {");
		group.add(fun0RadioButton);
		fun0RadioButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				integrationGraphics.odeNr = 0;
				updateGraphics();
			}
		});

		fun0RadioButton.setSelected(true);
		rightPanel.add(fun0RadioButton, c);
		c.gridy++;
		rightPanel.add(new JLabel("          vstart_x,"), c);
		c.gridy++;
		rightPanel.add(new JLabel("          vstart_y - g*t  }"), c);
		c.gridy++;

		JRadioButton fun1RadioButton = new JRadioButton("f(t,x) = {");
		group.add(fun1RadioButton);
		fun1RadioButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				integrationGraphics.odeNr = 1;
				updateGraphics();
			}
		});
		rightPanel.add(fun1RadioButton, c);
		c.gridy++;
		rightPanel.add(new JLabel("          vstart_x + p_x*p_x*0.005,"), c);
		c.gridy++;
		rightPanel.add(new JLabel(
				"          vstart_y - g*t -p_x*p_x*t*t*0.001 }"), c);
		c.gridy++;

		JRadioButton fun2RadioButton = new JRadioButton("f(t,x) = {");
		group.add(fun2RadioButton);
		fun2RadioButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				integrationGraphics.odeNr = 2;
				updateGraphics();
			}
		});
		rightPanel.add(fun2RadioButton, c);
		c.gridy++;
		rightPanel.add(new JLabel("          vstart_x + t*4.0,"), c);
		c.gridy++;
		rightPanel.add(new JLabel("          vstart_y + 0.1*t^4 - 8.0*t^2 }"),
				c);
		c.gridy++;

		JRadioButton fun3RadioButton = new JRadioButton("f(t,x) = {");
		group.add(fun3RadioButton);
		fun3RadioButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				integrationGraphics.odeNr = 3;
				updateGraphics();
			}
		});

		rightPanel.add(fun3RadioButton, c);
		c.gridy++;
		rightPanel.add(new JLabel("          p_y-1,"), c);
		c.gridy++;
		rightPanel.add(new JLabel(
				"          -100*(pos_x-10) - 101 * pos_y + 1}"), c);
		c.gridy++;

		JPanel somePanel = new JPanel();
		somePanel.setLayout(new GridLayout(6, 2));

		/*
		 * TIMESTEP
		 */
		spinnerTimestepSize = new JSpinner();

		modelSpinnerTimestepSize = new SpinnerNumberModel(0.5, // initial value
				0.01, // min
				1000.0, // max
				0.01 // step
		);

		integrationGraphics.timestepSize = modelSpinnerTimestepSize.getNumber()
				.doubleValue();
		spinnerTimestepSize.setModel(modelSpinnerTimestepSize);

		spinnerTimestepSize.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				integrationGraphics.timestepSize = modelSpinnerTimestepSize
						.getNumber().doubleValue();
				updateGraphics();
			}
		});
		somePanel.add(spinnerTimestepSize);
		somePanel.add(new JLabel("Timestep Size"));

		/*
		 * VELOCITY X
		 */
		spinnerVelocityX = new JSpinner();

		modelSpinnerVelocityX = new SpinnerNumberModel(30.0, // initial value
				0.0, // min
				1000.0, // max
				0.5 // step
		);

		integrationGraphics.startVelocity[0] = modelSpinnerVelocityX
				.getNumber().doubleValue();
		spinnerVelocityX.setModel(modelSpinnerVelocityX);

		spinnerVelocityX.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				integrationGraphics.startVelocity[0] = modelSpinnerVelocityX
						.getNumber().doubleValue();
				updateGraphics();
			}
		});
		somePanel.add(spinnerVelocityX);
		somePanel.add(new JLabel("Initial Velocity X"));

		/*
		 * VELOCITY Y
		 */
		spinnerVelocityY = new JSpinner();

		modelSpinnerVelocityY = new SpinnerNumberModel(80.0, // initial value
				0.0, // min
				1000.0, // maY
				0.5 // step
		);

		integrationGraphics.startVelocity[1] = modelSpinnerVelocityY
				.getNumber().doubleValue();
		spinnerVelocityY.setModel(modelSpinnerVelocityY);

		spinnerVelocityY.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				integrationGraphics.startVelocity[1] = modelSpinnerVelocityY
						.getNumber().doubleValue();
				updateGraphics();
			}
		});
		somePanel.add(spinnerVelocityY);
		somePanel.add(new JLabel("Initial Velocity Y"));

		/*
		 * StartPosition X
		 */
		spinnerStartPositionX = new JSpinner();

		modelSpinnerStartPositionX = new SpinnerNumberModel(-350.0, // initial
																	// value
				-400.0, // min
				400.0, // max
				1 // step
		);

		integrationGraphics.startPosition[0] = modelSpinnerStartPositionX
				.getNumber().doubleValue();
		spinnerStartPositionX.setModel(modelSpinnerStartPositionX);

		spinnerStartPositionX.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				integrationGraphics.startPosition[0] = modelSpinnerStartPositionX
						.getNumber().doubleValue();
				updateGraphics();
			}
		});
		somePanel.add(spinnerStartPositionX);
		somePanel.add(new JLabel("Initial StartPosition X"));

		/*
		 * StartPosition Y
		 */
		spinnerStartPositionY = new JSpinner();

		modelSpinnerStartPositionY = new SpinnerNumberModel(0.0, // initial
																	// value
				-50.0, // min
				550.0, // maY
				1 // step
		);

		integrationGraphics.startPosition[1] = modelSpinnerStartPositionY
				.getNumber().doubleValue();
		spinnerStartPositionY.setModel(modelSpinnerStartPositionY);

		spinnerStartPositionY.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				integrationGraphics.startPosition[1] = modelSpinnerStartPositionY
						.getNumber().doubleValue();
				updateGraphics();
			}
		});
		somePanel.add(spinnerStartPositionY);
		somePanel.add(new JLabel("Initial StartPosition Y"));

		rightPanel.add(somePanel, c);
		c.gridy++;

		rightPanel.add(new JLabel(""), c);
		c.gridy++;

		rightPanel.add(new JSeparator(), c);
		c.gridy++;

		rightPanel.add(new JLabel(""), c);
		c.gridy++;

		/*
		 * EULER
		 */
		checkBoxEuler = new JCheckBox("Euler (red)", true);
		integrationGraphics.eulerIntegration = true;
		checkBoxEuler.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent arg0) {
				integrationGraphics.eulerIntegration = checkBoxEuler
						.isSelected();
				updateGraphics();
			}
		});
		c.gridy++;
		rightPanel.add(checkBoxEuler, c);

		textAreaEuler = new JTextArea();
		textAreaEuler.setEnabled(false);
		c.gridy++;
		rightPanel.add(textAreaEuler, c);

		/*
		 * HEUN
		 */
		checkBoxHeun = new JCheckBox("Heun (Blue)", true);
		integrationGraphics.heunIntegration = true;
		checkBoxHeun.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				integrationGraphics.heunIntegration = checkBoxHeun.isSelected();
				updateGraphics();
			}
		});
		c.gridy++;
		rightPanel.add(checkBoxHeun, c);

		textAreaHeun = new JTextArea();
		textAreaHeun.setEnabled(false);
		c.gridy++;
		rightPanel.add(textAreaHeun, c);

		/*
		 * RUNGE KUTTA
		 */
		checkBoxRungeKutta = new JCheckBox("Runge Kutta (Green)", true);
		checkBoxRungeKutta.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent arg0) {
				integrationGraphics.rungeKuttaIntegration = checkBoxRungeKutta
						.isSelected();
				updateGraphics();
			}
		});
		c.gridy++;
		rightPanel.add(checkBoxRungeKutta, c);

		textAreaRungeKutta = new JTextArea();
		textAreaRungeKutta.setEnabled(false);
		c.gridy++;
		rightPanel.add(textAreaRungeKutta, c);

		setLayout(new BorderLayout());
		add(rightPanel, BorderLayout.EAST);
		add((JPanel) integrationGraphics, BorderLayout.WEST);

		updateGraphics();
	}
}
