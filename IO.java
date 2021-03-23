package ABB;

import java.io.*;
import java.util.StringTokenizer;

public class IO 
{
    public static long[][] readMatrix(String path)
    {
        File fichero = new File(path);
        StringTokenizer token;
        String line, s;
        int filas=0, columnas=0, fila=0, col=0, max=0;
        long num;
        long[][] matriz = null;     
        
        try
        {
            FileReader fr = new FileReader(fichero);
            BufferedReader br = new BufferedReader(fr);
            
            line = br.readLine();
            token = new StringTokenizer(line);
            filas = Integer.parseInt(token.nextToken());
            columnas = Integer.parseInt(token.nextToken());
            matriz = new long[filas][columnas];
            
            while( (line = br.readLine()) != null)
            {
                token = new StringTokenizer(line," ,{}");
                
                while( token.hasMoreElements() )
                {
                    s = token.nextToken();
                    if( s.length() > max)
                        max = s.length();
                    
                    num = Long.parseLong( s );
                    matriz[fila][col] = num;
                    
                    col = (++col)%columnas;
                    if(col == 0)
                        fila++;
                }
            }            
        }catch(IOException ioe)
        {
            System.err.println( ioe.getMessage() );
        }
        
        return matriz;
    }
    
    public static void writeJPP( PublicParam pp, String folderpath, String instanceName )
    {
        try
        {
            String file = folderpath + "\\PublicParameters_" + instanceName + ".jpp";
            FileOutputStream fos = new FileOutputStream( file );
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            
            oos.writeObject( pp );
            oos.close();
            fos.close();
            String out = "\n >> file 'PublicParameters_" + instanceName + ".jpp' generated";

        }catch(IOException ioe){}
    }
    
    public static PublicParam readJPP( String path)
    {
        PublicParam pp = null;
        
        try
        {
            FileInputStream fis = new FileInputStream( path );
            ObjectInputStream ois = new ObjectInputStream( fis );
            
            pp = (PublicParam)ois.readObject();
            ois.close();
            fis.close();
        }catch(IOException ioe){}
        catch(ClassNotFoundException cnfe)
        {
            System.out.println("Error : " + cnfe.getMessage() );
        }
        
        return pp;
    }
    
    public static void writePrivateKey( PrivateKey prvkey, String pathToFile )
    {
        try
        {
            FileOutputStream fos = new FileOutputStream( pathToFile );
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            
            System.out.println(" writing object");
            oos.writeObject( prvkey );
            oos.close();
            fos.close();
            System.out.println(" prvkey done ");
        }catch(IOException ioe){}

    }
    
    public static PrivateKey readPrivateKey( String path)
    {
        PrivateKey prvkey = null;
        
        try
        {
            FileInputStream fis = new FileInputStream( path );
            ObjectInputStream ois = new ObjectInputStream( fis );
            
            prvkey = (PrivateKey)ois.readObject();
            ois.close();
            fis.close();
        }catch(IOException ioe){}
        catch(ClassNotFoundException cnfe)
        {
            System.out.println("Error : " + cnfe.getMessage() );
        }
        
        return prvkey;
    }
    
    public static void saveGSO(String pathToFile, realVector[] gso)
    {
        ObjectOutputStream output = null;
        
        try
        {
            output = new ObjectOutputStream( new FileOutputStream( pathToFile ) );
        }
        catch( IOException ioe)
        {
            System.err.println(" Error opening file : " + pathToFile);
        }
        
        int i;
        
        for(i=0; i<gso.length; i++)
        {
            try
            {
                output.writeObject( gso[i] );
            }
            catch( IOException ioe )
            {
                System.err.println(" Error writing to file : " + pathToFile );
            }
        }
        
        try
        {
            if( output != null )
                output.close();
        }
        catch( IOException ioe)
        {
            System.err.println(" Error closing file : " + pathToFile );
            System.exit( 1 );
        }
    }
    
    public static realVector[] readGSO(String pathToFile, int m)
    {
        ObjectInputStream input = null;
        realVector[] gso = null;
        int i=0;
        
        try
        {
            input = new ObjectInputStream( new FileInputStream( pathToFile ));
            gso = new realVector[m];
        }
        catch(IOException ioe)
        {
            System.err.println(" Error opening file : " + pathToFile );
        }
        
        try
        {
            while( true )
            {
                gso[i++] = (realVector)input.readObject();
            }
        }
        catch( EOFException eofe)
        {
            
        }
        catch( ClassNotFoundException cnfe)
        {
            System.err.println(" Unable to create object from file : " + pathToFile);
        }
        catch( IOException ioe)
        {
            System.err.println(" Error during read from file : " + pathToFile);
        }
        
        try
        {
            if( input != null )
                input.close();
        }
        catch( IOException ioe )
        {
            System.err.println(" Error closing file : " + pathToFile);
        }

        return gso;
    }
    
    
}
