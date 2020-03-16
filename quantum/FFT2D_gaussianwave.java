import org.opensourcephysics.controls.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.numerics.*;

public class FFT2D_gaussianwave extends AbstractCalculation {
  FFT2DFrame kFrame = new FFT2DFrame("k_x", "k_y", "2D FFT");
  Complex2DFrame xFrame = new Complex2DFrame("x","y","Complex Field");

  public void calculate() {
    xFrame.clearData();
    String fStr = control.getString("f(|r|)");
    Function f = null;
    try{
      f = new ParsedFunction(fStr, "r");
    } catch(ParserException ex){
      control.println("Error parsing function string: "+fStr);
      return;
    }
    int xMode = control.getInt("x mode"),
        yMode = control.getInt("y mode");
    double xmin = control.getDouble("xmin");
    double xmax = control.getDouble("xmax");
    int nx = control.getInt("Nx");
    double ymin = control.getDouble("ymin");
    double ymax = control.getDouble("ymax");
    int ny = control.getInt("Ny");

    // data stored in row-major format
    double[] zdata = new double[2*nx*ny];
    xFrame.setPreferredMinMax(xmin,xmax,ymin,ymax);
    xFrame.resizeGrid(nx,ny);
    xFrame.setAll(zdata);
    for(int iy=0; iy<ny; iy++){
      // offset to beginning of a row; each row length = nx
      int offset = 2*iy*nx;
      double y = xFrame.indexToY(iy);
      for(int ix=0; ix<nx; ix++){
        double x = xFrame.indexToX(ix);
        // f(z) = e^(i*xmode*x)e^(i*ymode*y)
        zdata[offset+2*ix] = // real part
           f.evaluate(x*x+y*y)*(Math.cos(xMode*x)*Math.cos(yMode*y)
          -Math.sin(xMode*x)*Math.sin(yMode*y));
        zdata[offset+2*ix+1] = // imag part
           f.evaluate(x*x+y*y)*(Math.sin(xMode*x)*Math.cos(yMode*y)
          +Math.cos(xMode*x)*Math.sin(yMode*y));
      }

    }
    xFrame.setAll(zdata);
    xFrame.setVisible(true);
    xFrame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
    kFrame.doFFT(zdata, nx, xmin, xmax, ymin, ymax);
  }

  public void reset() {
    control.setValue("f(|r|)", "exp(-r*r)");
    control.setValue("x mode",0);
    control.setValue("y mode", 5);
    control.setValue("xmin", "-pi");
    control.setValue("xmax", "pi");
    control.setValue("ymin", "-pi");
    control.setValue("ymax", "pi");
    control.setValue("Nx",32);
    control.setValue("Ny",32);
  }

  public static void main(String[] args){
    CalculationControl.createApp(new FFT2D_gaussianwave());
  }
}
