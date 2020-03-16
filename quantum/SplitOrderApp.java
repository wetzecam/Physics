import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.numerics.*;

public class SplitOrderApp extends AbstractSimulation {
  Complex2DFrame xFrame = new Complex2DFrame("x", "y", "Wave function");
  //Complex2DFrame pFrame = new Complex2DFrame("x", "y", "Wave function");
  SplitOrder wavefunction;
  double time;

  public SplitOrderApp(){
    xFrame.setShowGrid(false);
    xFrame.setSize(600,600);
    return;
  }

  public void initialize() {
    time =0;
    xFrame.setMessage("t = "+0);

    String fStr = control.getString("f(|r|)");
    Function f = null;
    try{
      f = new ParsedFunction(fStr, "r");
    } catch(ParserException ex){
      control.println("Error parsing function string: "+fStr);
      return;
    }
    int xMode = control.getInt("x mode");
    int yMode = control.getInt("y mode");
    double xmin = control.getDouble("xmin");
    double xmax = control.getDouble("xmax");
    double ymin = control.getDouble("ymin");
    double ymax = control.getDouble("ymax");
    int nx = control.getInt("Nx");
    int ny = control.getInt("Ny");
    double xinit = control.getDouble("x init");
    double yinit = control.getDouble("y init");


    wavefunction = new SplitOrder(f, xMode, yMode, xmin, xmax, ymin, ymax, nx, ny, xinit, yinit);

    xFrame.clearData();
    xFrame.setPreferredMinMax(xmin,xmax,ymin,ymax);
    xFrame.resizeGrid(nx,ny);
    xFrame.setAll(wavefunction.Psi);


  }

  public void doStep(){
    time += wavefunction.step();
    time += wavefunction.step();
    xFrame.setAll(wavefunction.Psi);
    xFrame.setVisible(true);
    xFrame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
    xFrame.setMessage("t = "+decimalFormat.format(time));
  }

  public void reset(){
    control.setValue("f(|r|)", "exp(-r*r/2)");
    control.setValue("x mode", 2);
    control.setValue("y mode", 5);
    control.setValue("xmin", "-4*pi");
    control.setValue("xmax",  "4*pi");
    control.setValue("ymin", "-4*pi");
    control.setValue("ymax",  "4*pi");
    control.setValue("Nx", 64);
    control.setValue("Ny", 64);
    control.setValue("x init", -5.0);
    control.setValue("y init", -3.0);

    // multiple computation steps per animatoin step
    setStepsPerDisplay(20);
    enableStepsPerDisplay(true);
    initialize();
  }

  public static void main(String[] args){
    SimulationControl.createApp(new SplitOrderApp());
  }
}
