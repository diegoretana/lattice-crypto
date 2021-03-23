package ABB;

import com.wolfram.jlink.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Date;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

import java.util.StringTokenizer;

public class AlwenPeikert
{
    private int m1;
    private int m2;
    
    private int l, r;
    private long q;
    private int i, j, k; // counters
    
    private LMatrix A;
    private LMatrix S;
    private LMatrix Hp;
    
    long[][] H;
    
    public AlwenPeikert( LMatrix A1, LMatrix orto, int r, int l, long q )
    {
        KernelLink ml = SystemParameters.ml;
        this.r = r;
        this.q = q;
        this.l = l;
        this.m1 = A1.columns;
        this.m2 = m1 * l;
        
        String cad;
        try
        {
            cad = "A = " + orto.transpose().getString();
            ml.evaluate(cad);
            ml.discardAnswer();
            cad = "{u,H} = HermiteDecomposition[A];";
            ml.evaluate(cad);
            ml.discardAnswer();
            cad = "H";
            ml.evaluate(cad);
            ml.waitForAnswer();
            H = ml.getLongArray2();
            
            for(i=0; i<H.length; i++)
                H[i][i] -= 1; // H'
        
            Hp = new LMatrix(H);
            Hp = Hp.transpose();

            H = Hp.M.clone();
        }
        catch(MathLinkException mle)
        {
            System.out.println(" Error : " + mle.getMessage() );
        
        }
    }
    
    public LMatrix setup( LMatrix A1, String folder, String instanceName )
    {
        return Algorithm1(A1, folder, instanceName );
    }
    
    
    private LMatrix Algorithm1( LMatrix A1, String folder, String instanceName )
    {
        LMatrix A2, U, G, R, P;
        
        U = Define_U();
        P = Define_P();
        R = Define_R();
        G = Define_G();
        
        A1.mult_const(-1);
        A2 = A1.ModMult(R.addition(G), q);
        
        A1.mult_const(-1); // reverse prior multiplication
        A = A1.augment(A2);
        
        long[][] c = new long[Hp.rows][Hp.columns];
        for(i=0; i<Hp.rows; i++)
            c[i][i] = 1;
        LMatrix C = new LMatrix(c);
        C.setIdentity();
        C.mult_const(-1);
        LMatrix M01 = ( R.mult(P) ).addition( C );        
        
        LMatrix M00 = ( G.addition(R) ).mult(U);
        LMatrix top  = M00.augment(M01);
        LMatrix down = U.augment(P);
        
        S = top.augmentDown(down);

        String file = folder + "\\Ta_" + instanceName + ".tdr";
        S.WriteLMatrix(file);
        String out = "\n >> file 'Ta_" + instanceName + ".tdr' generated ";
        
        return A;
    }
    
    public LMatrix Define_P()
    {
        long[][] P = new long[m2][m1];
        
        for(j=0; j<m1; j++)
            P[(j+1)*l-1][j] = 1;
        
        return new LMatrix(P);
    }
    
    public LMatrix Define_R()
    {
        SecureRandom sr = null;
        
        try
        {
            sr = SecureRandom.getInstance("SHA1PRNG");
            sr.setSeed(new Date().getTime());
        }catch(NoSuchAlgorithmException nsae){}
        
        long[][] R = new long[m1][m2];
        
        int val;
        for(i=0; i<m1; i++)
            for(j=0; j<m2; j++)
            {
                val = sr.nextInt(4);
                
                try{
                switch(val)
                {
                    case 0 : R[i][j] = 1; break;
                    case 1 : 
                    case 2 : R[i][j] = 0; break;
                    case 3 : R[i][j] = -1; 
                }
                }catch(ArrayIndexOutOfBoundsException e){}
            }
        
        return new LMatrix(R);
    }
    
    public LMatrix Define_U()
    {
        long[][] U = new long[m2][m2];
        
        for(i=0; i<m2; i++)
            U[i][i] = 1;
        
        for(k=0; k<m1; k++)
            for(i=0; i<= l-2; i++)
                U[k*l+i][k*l+i+1] = -r;

        return new LMatrix(U);
    }

    private LMatrix Define_G()
    {
        long[][] G = new long[m1][m2];
        
        long h;
        for(i=1; i<=m1; i++)
            for(j=0; j<m1; j++)
            {
                h = H[j][i-1];
                G[j][i*l-1] = h;
                
                for(k=1; k<l; k++)
                    G[j][i*l-1-k] = (long)Math.floor( h/Math.pow(r, k) );
            }
        
        return new LMatrix(G);
    }
    
    public LMatrix getA()
    {
        return A;
    }
    
    public LMatrix getS()
    {
        return S;
    }
    
    public static LMatrix readLMatrix(String path)
    {
        File fichero = new File(path);
        StringTokenizer token;
        String line, s;
        int filas=0, columnas=0, fila=0, col=0, max=0;
        long[][] matriz = null;
        long num;
        
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
        
        
        max++;
        for(fila=0; fila<filas; fila++)
        {
            //System.out.print((fila+1) + "  ");
            for(col=0; col<columnas; col++)
            {
                s = String.valueOf(matriz[fila][col]);
                while( s.length() < max )
                    s = " " + s;
                
                System.out.print(s);
            }
            System.out.println("");
        }
       
        
        return new LMatrix(matriz);
    }
}
