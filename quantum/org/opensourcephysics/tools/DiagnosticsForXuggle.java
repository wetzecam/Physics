package org.opensourcephysics.tools;

import java.io.File;
import java.lang.reflect.Method;
import java.text.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.logging.Level;

import javax.swing.*;

import org.opensourcephysics.controls.OSPLog;
import org.opensourcephysics.controls.XML;
import org.opensourcephysics.display.OSPRuntime;
import org.opensourcephysics.media.core.VideoIO;

/**
 * Checks to see if Xuggle is installed and working.
 * 
 * @author Wolfgang Christian
 * @author Douglas Brown
 * @version 1.0
 */
public class DiagnosticsForXuggle {
	
  @SuppressWarnings("javadoc")
	public static final String XUGGLE_URL ="http://www.compadre.org/osp/items/detail.cfm?ID=11606"; //$NON-NLS-1$
	static String newline = System.getProperty("line.separator", "\n"); //$NON-NLS-1$ //$NON-NLS-2$
	static String[] xuggleJarNames = new String[] {"xuggle-xuggler.jar", "logback-core.jar",  //$NON-NLS-1$ //$NON-NLS-2$
			"logback-classic.jar", "slf4j-api.jar"}; //$NON-NLS-1$ //$NON-NLS-2$
	static boolean is64BitVM;
	static String xuggleHome, javaExtDirectory;
	static String pathEnvironment, pathValue;
	static File[] xuggleHomeJars, javaExtensionJars;
	
	static{  // added by W. Christian
    try {
      String name = "org.opensourcephysics.media.xuggle.XuggleIO"; //$NON-NLS-1$
      Class<?> xuggleClass = Class.forName(name);
      Method method=xuggleClass.getMethod("registerWithVideoIO"); //$NON-NLS-1$
      method.invoke(null, (Object[])null);
		} catch(Exception ex) {
		} catch(Error err) {
		}
    
    String s = System.getProperty("java.vm.name"); //$NON-NLS-1$
  	s += "-"+System.getProperty("os.arch"); //$NON-NLS-1$ //$NON-NLS-2$
  	s += "-"+System.getProperty("sun.arch.data.model"); //$NON-NLS-1$ //$NON-NLS-2$
  	is64BitVM = s.indexOf("64")>-1; //$NON-NLS-1$
  	
		// get XUGGLE_HOME
    xuggleHome = System.getenv("XUGGLE_HOME"); //$NON-NLS-1$	
    xuggleHomeJars = new File[xuggleJarNames.length];
    javaExtensionJars = new File[xuggleJarNames.length];
  }

	private DiagnosticsForXuggle() {
	}
	
	/**
	 * Displays the About Xuggle dialog. If working correctly, shows version, etc.
	 * If not working, shows a diagnostic message.
	 */
	public static void aboutXuggle() {

		int status = getStatusCode();
		
		if (OSPLog.getLevelValue()<=Level.CONFIG.intValue()) {
			// log XUGGLE_HOME and PATH environment variables
			String xuggleHome = System.getenv("XUGGLE_HOME"); //$NON-NLS-1$
		  OSPLog.config("XUGGLE_HOME = "+xuggleHome); //$NON-NLS-1$
		  OSPLog.config(pathEnvironment+" = "+pathValue); //$NON-NLS-1$
			
			// log xuggle home jars
			File[] xuggleJars = getXuggleJarFiles(xuggleHome+"/share/java/jars"); //$NON-NLS-1$
			boolean hasAllHomeJars = xuggleJars[0]!=null;
			for (int i=1; i< xuggleJars.length; i++) {
				hasAllHomeJars = hasAllHomeJars && xuggleJars[i]!=null;
			}
			SimpleDateFormat sdf = new SimpleDateFormat();
			String[] jarDates = new String[xuggleJarNames.length];
			for (int i=0; i< jarDates.length; i++) {
				jarDates[i] = xuggleJars[i]==null? "": " modified "+sdf.format(xuggleJars[i].lastModified()); //$NON-NLS-1$ //$NON-NLS-2$
			}
			String[] jarSizes = new String[xuggleJarNames.length];
			for (int i=0; i< jarSizes.length; i++) {
				jarSizes[i] = xuggleJars[i]==null? "": " (file size "+(xuggleJars[i].length()/1024)+"kB) "; //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
			}
			String fileData = "XUGGLE_HOME Xuggle files: ";			 //$NON-NLS-1$
			for (int i=0; i< jarSizes.length; i++) {
				if (i>0) fileData+=", "; //$NON-NLS-1$
				fileData+=xuggleJarNames[i]+" "+jarSizes[i]+xuggleJars[i]+jarDates[i]; //$NON-NLS-1$
			}
		  OSPLog.config(fileData);
		    
	    // log current java VM and extension jars
		  String javaHome = System.getProperty("java.home");	//$NON-NLS-1$
		  String bitness = is64BitVM? "(64-bit): ": "(32-bit): "; //$NON-NLS-1$ //$NON-NLS-2$
		  OSPLog.config("Java VM "+bitness+javaHome); //$NON-NLS-1$
			xuggleJars = getJavaExtensionJars();
			for (int i=0; i< jarDates.length; i++) {
				jarDates[i] = xuggleJars[i]==null? "": " modified "+sdf.format(xuggleJars[i].lastModified()); //$NON-NLS-1$ //$NON-NLS-2$
			}
			for (int i=0; i< jarSizes.length; i++) {
				jarSizes[i] = xuggleJars[i]==null? "": " (file size "+(xuggleJars[i].length()/1024)+"kB) "; //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
			}
			fileData = "Java extension Xuggle files: ";			 //$NON-NLS-1$
			for (int i=0; i< jarSizes.length; i++) {
				if (i>0) fileData+=", "; //$NON-NLS-1$
				fileData+=xuggleJarNames[i]+" "+jarSizes[i]+xuggleJars[i]+jarDates[i]; //$NON-NLS-1$
			}
		  OSPLog.config(fileData);
		  
	    // log tracker home jars on Windows
		  if (OSPRuntime.isWindows()) {
				String trackerHome = System.getenv("TRACKER_HOME"); //$NON-NLS-1$
			  OSPLog.config("TRACKER_HOME = "+trackerHome); //$NON-NLS-1$
			  xuggleJars =  getXuggleJarFiles(trackerHome);
				for (int i=0; i< jarDates.length; i++) {
					jarDates[i] = xuggleJars[i]==null? "": " modified "+sdf.format(xuggleJars[i].lastModified()); //$NON-NLS-1$ //$NON-NLS-2$
				}
				for (int i=0; i< jarSizes.length; i++) {
					jarSizes[i] = xuggleJars[i]==null? "": " (file size "+(xuggleJars[i].length()/1024)+"kB) "; //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
				}
				fileData = "TRACKER_HOME Xuggle files: ";			 //$NON-NLS-1$
				for (int i=0; i< jarSizes.length; i++) {
					if (i>0) fileData+=", "; //$NON-NLS-1$
					fileData+=xuggleJarNames[i]+" "+jarSizes[i]+xuggleJars[i]+jarDates[i]; //$NON-NLS-1$
				}
			  OSPLog.config(fileData);
		  }

		}
			
		// display appropriate dialog		
		if (status==0) { // xuggle working correctly
			String unknown = XuggleRes.getString("Xuggle.Dialog.Unknown"); //$NON-NLS-1$
			String version = unknown;
			try {
				String name = "com.xuggle.xuggler.Version"; //$NON-NLS-1$
		    Class<?> xuggleClass = Class.forName(name);
		    Method method = xuggleClass.getMethod("getVersionString"); //$NON-NLS-1$
		    version=(String)method.invoke(null, (Object[])null);
			} catch(Exception ex) {
			} catch(Error err) {
	    }			
			xuggleHome = " " + xuggleHome + newline; //$NON-NLS-1$
			String fileInfo = newline;
			String path = " "+unknown; //$NON-NLS-1$
			if (javaExtensionJars[0]!=null) {
				DateFormat format = DateFormat.getDateInstance(DateFormat.SHORT);
				Date date = new Date(javaExtensionJars[0].lastModified());
				fileInfo = " (" + format.format(date) + ")" + newline; //$NON-NLS-1$ //$NON-NLS-2$
				path = " " + javaExtensionJars[0].getAbsolutePath(); //$NON-NLS-1$
			}
			JOptionPane.showMessageDialog(null,
					XuggleRes.getString("Xuggle.Dialog.AboutXuggle.Message.Version") //$NON-NLS-1$
							+ " " + version + fileInfo //$NON-NLS-1$
							+ XuggleRes.getString("Xuggle.Dialog.AboutXuggle.Message.Home") //$NON-NLS-1$
							+ xuggleHome 
							+ XuggleRes.getString("Xuggle.Dialog.AboutXuggle.Message.Path") //$NON-NLS-1$ 
							+ path,
					XuggleRes.getString("Xuggle.Dialog.AboutXuggle.Title"), //$NON-NLS-1$
					JOptionPane.INFORMATION_MESSAGE);
		}
		
		else { // xuggle not working
    	String[] diagnostic = getDiagnosticMessage(status);
    	Box box = Box.createVerticalBox();
    	box.setBorder(BorderFactory.createEmptyBorder(0, 0, 10, 0));
    	for (String line: diagnostic) {    			
  			box.add(new JLabel(line));
  		}      	
    	JOptionPane.showMessageDialog(null, box,
    			XuggleRes.getString("Xuggle.Dialog.BadXuggle.Title"),  //$NON-NLS-1$
    			JOptionPane.WARNING_MESSAGE);
		}

	}
	
	/**
	 * Gets the xuggle jar files (named in xuggleJarNames) found in a given directory.
	 * Always returns the array, but individual elements may be null.
	 * @param dir the directory
	 * @return the array of jar files found
	 */
	public static File[] getXuggleJarFiles(String dir) {
		// look for xuggle-xuggler and support jars in the directory
		File[] jarFiles = new File[xuggleJarNames.length];
		for (int i=0; i< jarFiles.length; i++) {
			String next = xuggleJarNames[i];
			File file = new File(dir, next);
			jarFiles[i] = file.exists()? file: null;
		}
		return jarFiles;
	}
	
	
	/**
	 * Gets the xuggle jar files (named in xuggleJarNames) found in the current Java extension directory.
	 * Always returns the array, but individual elements may be null.
	 * @return the array of jar files found
	 */
	public static File[] getJavaExtensionJars() {
				
		File[] xuggleFiles = new File[xuggleJarNames.length];
		
		// look for xuggle jars in system extensions directories
		String extFolders = XML.forwardSlash(System.getProperty("java.ext.dirs")); //$NON-NLS-1$
		String separator = System.getProperty("path.separator"); //$NON-NLS-1$
		int n = extFolders.indexOf(separator);
		if (n==-1) { // no separators, so single path
			javaExtDirectory = extFolders;
			xuggleFiles = getXuggleJarFiles(extFolders);
			if (xuggleFiles[0]!=null) {
				return xuggleFiles;
			}
		}
		else {
			String dir = extFolders;
			javaExtDirectory = null;
			while (xuggleFiles[0]==null && n > -1) {
				dir = extFolders.substring(0, n);
				if (javaExtDirectory==null)
					javaExtDirectory = dir;	// first one in list by default		
				extFolders = extFolders.substring(n + 1);
				xuggleFiles = getXuggleJarFiles(dir);
				if (xuggleFiles[0]!=null) {
					javaExtDirectory = dir;
					return xuggleFiles;
				}
				n = extFolders.indexOf(separator);
			}			
		}
		
		return xuggleFiles;
	}
	
	
	/**
	 * Gets a status code that identifies the current state of the Xuggle video engine.
	 * Codes are:
	 * 		0	working correctly
	 * 		1 not installed (XUGGLE_HOME==null, incomplete xuggle jars in current java extensions directory)
	 * 		2	needs reboot/login? (XUGGLE_HOME==null, complete xuggle jars in current java extensions directory)
	 * 		3 XUGGLE_HOME defined but incomplete xuggle jars in xuggle home
	 * 		4 XUGGLE_HOME complete, but incorrect "PATH", "DYLD_LIBRARY_PATH", or "LD_LIBRARY_PATH"
	 * 		5 XUGGLE_HOME complete, but missing xuggle jars in current java extensions directory
	 * 		6 XUGGLE_HOME complete, but mismatched xuggle-xuggler.jar in current java extensions directory
	 * 		7 XUGGLE_HOME complete, but wrong Java VM bitness
	 * 	 -1 none of the above
	 * 
	 * @return status code
	 */
	public static int getStatusCode() {
		javaExtensionJars = getJavaExtensionJars();		
		pathEnvironment = OSPRuntime.isWindows()? "PATH": OSPRuntime.isMac()? "DYLD_LIBRARY_PATH": "LD_LIBRARY_PATH"; //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
		pathValue = System.getenv(pathEnvironment);
		
		// return 0 if working correctly
		if (VideoIO.getVideoType("xuggle", null)!=null) return 0; //$NON-NLS-1$
		
		boolean completeExt = javaExtensionJars[0]!=null;
		for (int i=1; i< javaExtensionJars.length; i++) {
			completeExt = completeExt && javaExtensionJars[i]!=null;
		}
		
		if (xuggleHome==null) {
			return completeExt? 2: 1;			
		}
		
		// get xuggle home jars
		xuggleHomeJars = getXuggleJarFiles(xuggleHome+"/share/java/jars"); //$NON-NLS-1$
		boolean completeHome = xuggleHomeJars[0]!=null;
		for (int i=1; i< xuggleHomeJars.length; i++) {
			completeHome = completeHome && xuggleHomeJars[i]!=null;
		}
		
		if (!completeHome) return 3;
		if (javaExtensionJars[0]==null) return 5;
		
		long homeLength = xuggleHomeJars[0].length();
		boolean mismatched = homeLength!=javaExtensionJars[0].length();
		if (mismatched) return 6;
		if (!completeExt) return 5;
		
		if (homeLength<1000000) { // older Xuggle version (probably 3.4) installed
			String folder = OSPRuntime.isWindows()? "/bin": "/lib"; //$NON-NLS-1$ //$NON-NLS-2$
			String xuggleLib = XML.forwardSlash(xuggleHome+folder);
			if (XML.forwardSlash(pathValue).indexOf(xuggleLib)==-1) return 4;			
			if (is64BitVM && OSPRuntime.isWindows()) return 7;
		}
		
		if (!is64BitVM && OSPRuntime.isMac()) return 7;
		
		return -1;
	}

	
	/**
	 * Gets a diagnostic message when Xuggle is not working.
	 * @param status the status code from getStatusCode() method
	 * @return an array strings containing the message lines
	 */
	public static String[] getDiagnosticMessage(int status) {
		
		if (status==0) return null;

		ArrayList<String> message = new ArrayList<String>();		
		switch(status) {
			case 1: // not installed (XUGGLE_HOME==null, incomplete xuggle jars in current java extensions directory)
				message.add(XuggleRes.getString("Xuggle.Dialog.NoXuggle.Message1")); //$NON-NLS-1$
	    	message.add(" "); //$NON-NLS-1$
	    	message.add(XuggleRes.getString("Xuggle.Dialog.AboutXuggle.Message.InstallerPath")); //$NON-NLS-1$
	    	message.add(XUGGLE_URL);
				break;
				
			case 2: // needs reboot/login? (XUGGLE_HOME==null, complete xuggle jars in current java extensions directory)
				message.add(XuggleRes.getString("Xuggle.Dialog.BadXuggle.Message")); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.NoXuggleHome.Message1")); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.NoXuggleHome.Message2")); //$NON-NLS-1$
	    	message.add(" "); //$NON-NLS-1$
	    	message.add(XuggleRes.getString("Xuggle.Dialog.AboutXuggle.Message.InstallerPath")); //$NON-NLS-1$
	    	message.add(XUGGLE_URL);
				break;
				
			case 3: // XUGGLE_HOME defined but incomplete xuggle jars in xuggle home
				message.add(XuggleRes.getString("Xuggle.Dialog.BadXuggle.Message")); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.IncompleteXuggle.Message1")); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.IncompleteXuggle.Message2")); //$NON-NLS-1$
	    	message.add(" "); //$NON-NLS-1$
	    	message.add(XuggleRes.getString("Xuggle.Dialog.AboutXuggle.Message.InstallerPath")); //$NON-NLS-1$
	    	message.add(XUGGLE_URL);
				break;
				
			case 4: // XUGGLE_HOME complete, but incorrect "PATH", "DYLD_LIBRARY_PATH", or "LD_LIBRARY_PATH"
				message.add(XuggleRes.getString("Xuggle.Dialog.MissingEnvironmentVariable.Message1")); //$NON-NLS-1$
				message.add("\""+pathEnvironment+"\" "+XuggleRes.getString("Xuggle.Dialog.MissingEnvironmentVariable.Message2")); //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
				message.add(XuggleRes.getString("Xuggle.Dialog.MissingEnvironmentVariable.Message3")); //$NON-NLS-1$
				break;
				
			case 5: // XUGGLE_HOME complete, but missing xuggle jars in current java extensions directory
				String missingJars = ""; //$NON-NLS-1$
				for (int i=0; i<xuggleJarNames.length; i++) {
					if (javaExtensionJars[i]==null) {
						if (missingJars.length()>1) missingJars += ", "; //$NON-NLS-1$
						missingJars += xuggleJarNames[i];
					}
				}
				String source = XML.forwardSlash(xuggleHome)+"/share/java/jars"; //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.MissingJarsInExt.Message1")); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.MissingJarsInExt.Message2")); //$NON-NLS-1$
	    	message.add(" "); //$NON-NLS-1$
	    	message.add(missingJars);
	    	message.add(" "); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.MissingJarsInExt.Message3")); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.MissingJarsInExt.Message4")); //$NON-NLS-1$
	    	message.add(" "); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.SourceDirectory.Message")+" "+source); //$NON-NLS-1$ //$NON-NLS-2$
				message.add(XuggleRes.getString("Xuggle.Dialog.ExtensionDirectory.Message")+" "+javaExtDirectory); //$NON-NLS-1$ //$NON-NLS-2$
				break;
				
			case 6: // XUGGLE_HOME complete, but mismatched xuggle-xuggler.jar in current java extensions directory
				source = XML.forwardSlash(xuggleHome)+"/share/java/jars"; //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.MismatchedJar.Message1")); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.MismatchedJar.Message2")); //$NON-NLS-1$
	    	message.add(" "); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.MismatchedJar.Message3")); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.MismatchedJar.Message4")); //$NON-NLS-1$
	    	message.add(" "); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.SourceDirectory.Message")+" "+source); //$NON-NLS-1$ //$NON-NLS-2$
				message.add(XuggleRes.getString("Xuggle.Dialog.ExtensionDirectory.Message")+" "+javaExtDirectory); //$NON-NLS-1$ //$NON-NLS-2$
				break;
				
			case 7: // XUGGLE_HOME complete, but wrong Java VM bitness
				if (OSPRuntime.isMac()) {
					// wrong VM on Mac
		    	message.add(XuggleRes.getString("Xuggle.Dialog.WrongVMMac.Message1")); //$NON-NLS-1$
		    	message.add(XuggleRes.getString("Xuggle.Dialog.WrongVMMac.Message2")); //$NON-NLS-1$
		    	message.add(" "); //$NON-NLS-1$
		    	message.add(XuggleRes.getString("Xuggle.Dialog.WrongVMMac.Message3")); //$NON-NLS-1$
		    	message.add(XuggleRes.getString("Xuggle.Dialog.WrongVMMac.Message4")); //$NON-NLS-1$
		    	message.add(XuggleRes.getString("Xuggle.Dialog.WrongVMMac.Message5")); //$NON-NLS-1$					
				}
				else {
					// wrong VM on Windows
		    	message.add(XuggleRes.getString("Xuggle.Dialog.WrongVMWindows.Message1")); //$NON-NLS-1$
		    	message.add(XuggleRes.getString("Xuggle.Dialog.WrongVMWindows.Message2")); //$NON-NLS-1$
		    	message.add(" "); //$NON-NLS-1$
		    	
					Collection<File> jreDirs = new ExtensionsManager().findPublicJREsForXuggle(is64BitVM);
					if (jreDirs.isEmpty()) {
			    	message.add(XuggleRes.getString("Xuggle.Dialog.NoVM.Message1")); //$NON-NLS-1$
			    	message.add(XuggleRes.getString("Xuggle.Dialog.NoVM.Message2")); //$NON-NLS-1$
			    	message.add(XuggleRes.getString("Xuggle.Dialog.NoVM.Message3")); //$NON-NLS-1$
			    	message.add(XuggleRes.getString("Xuggle.Dialog.NoVM.Message4")); //$NON-NLS-1$
					}
					else {
			    	message.add(XuggleRes.getString("Xuggle.Dialog.SetVM.Message1")); //$NON-NLS-1$
			    	message.add(XuggleRes.getString("Xuggle.Dialog.SetVM.Message2")); //$NON-NLS-1$
			    	message.add(" "); //$NON-NLS-1$
						for (File next: jreDirs) {
				    	message.add(next.getPath());
						}
					}
				}
				break;
				
			default: // none of the above
				message.add(XuggleRes.getString("Xuggle.Dialog.BadXuggle.Message")); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.UnknownProblem.Message")); //$NON-NLS-1$
				message.add(XuggleRes.getString("Xuggle.Dialog.NoXuggleHome.Message2")); //$NON-NLS-1$
	    	message.add(" "); //$NON-NLS-1$
	    	message.add(XuggleRes.getString("Xuggle.Dialog.AboutXuggle.Message.InstallerPath")); //$NON-NLS-1$
	    	message.add(XUGGLE_URL);
				
		}
		
		return message.toArray(new String[message.size()]);
	}
	
	/**
	 * Tests this class.
	 * @param args ignored
	 */
	public static void main(String[] args) {
		aboutXuggle();
  }
}



