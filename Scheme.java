package ABB;

import com.wolfram.jlink.*;
import java.util.StringTokenizer;
import java.math.BigDecimal;
import java.nio.ByteBuffer;

public class Scheme
{
    private LMatrix Ta;
    private LMatrix Fid = null;
    private LMatrix FidT = null;
    
    public static int n;
    public int m;
    public int q;
    
    public Scheme(PublicParam PP)
    {
        this.n = PP.n;
        this.q = PP.q;
        this.m = PP.m;
        SystemParameters.q = PP.q;
        SystemParameters.n = PP.n;
    }
    
    public Scheme(){}
    
    public void Setup(int n, String folder, String pathToKernel, String instanceName )
    {
        KernelLink ml = null;
        String out="";
        
        try 
        { // 'C:\\Program Files\\Wolfram Research\\Mathematica\\9.0\\MathKernel.exe'
            ml = MathLinkFactory.createKernelLink("-linkmode launch -linkname " + "'"+pathToKernel+"'");
        }catch(MathLinkException mle)
        {
            System.out.println("Fatal error opening link : " + mle.getMessage());
        }
        
        String expression;
        try
        {
            SystemParameters.ml = ml;
            this.n = n;
            SystemParameters.n = n;

            ml.discardAnswer();
            expression = "n = " + String.valueOf( n ) + ";";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "omega[n_] := Log2[n];";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "m = 6*n";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "q = NextPrime[m^(2.5)*omega[n]]";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "delta = 0;";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "While[n^delta <= Ceiling[Log2[q]],"
                    + "delta = delta + 0.001;"
                    + "m = Ceiling[6*n^(1+delta)];"
                    + "q = NextPrime[m^(2.5)*omega[n]]];";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "q";
            ml.evaluate(expression);
            ml.waitForAnswer();
            this.q = ml.getInteger();
            SystemParameters.q = this.q;
            out += "\n >> q : " + String.valueOf(this.q);
            
            expression = "Needs[\"FiniteFields`\"]";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "ToString[IrreduciblePolynomial[x,q,n], InputForm]";
            ml.evaluate(expression);
            ml.waitForAnswer();
            SystemParameters.irreduciblePoly = ml.getString();
            out += "\n >> polynomial : " + SystemParameters.irreduciblePoly;
            
            expression = "alpha = N[1/(m^2*omega[n])]";
            ml.evaluate(expression);
            ml.waitForAnswer();
            SystemParameters.alpha = ml.getDouble();
            out += "\n >> alpha : " + String.valueOf(SystemParameters.alpha);
            
            expression = "sigma = N[m*omega[n]]";
            ml.evaluate(expression);
            ml.waitForAnswer();
            SystemParameters.sigma = new BigDecimal( ml.getDouble() );
            out += "\n >> sigma : " + String.valueOf(SystemParameters.sigma);
            
            String file;
        
            LMatrix A0 = TrapGen( folder, instanceName );
            
            file = folder + "\\A0_" + instanceName + ".pp";
            A0.WriteLMatrix( file );
            out += "\n >> file 'A0_" + instanceName + ".pp' generated";

            LMatrix A1 = RandomGen.UniformRandom_Zq_nxm(n, m);
            file = folder + "\\A1_" + instanceName + ".pp";
            A1.WriteLMatrix( file );
            out += "\n >> file 'A1_" + instanceName + ".pp' generated";

            LMatrix B  = RandomGen.UniformRandom_Zq_nxm(n, m);
            file = folder + "\\B_" + instanceName + ".pp";
            B.WriteLMatrix( file );
            out += "\n >> file 'B_" + instanceName + ".pp' generated";

            PublicParam pp = new PublicParam(A0, A1, B, n, m, q, SystemParameters.irreduciblePoly, SystemParameters.alpha, SystemParameters.sigma );
            IO.writeJPP( pp, folder, instanceName );
            
        }catch(MathLinkException mle)
        {
            System.err.println(" Error : MathLinkException computing Setup " + 
                    mle.getMessage() );
            if(ml != null)
                ml.close();
        }

        ml.close();
        System.out.println(" process finished.");
    }

    public LMatrix TrapGen( String folder, String instanceName )
    {
        KernelLink ml = SystemParameters.ml;
        int l, r = SystemParameters.apr;
        LMatrix A = null;
        
        String expression;
        try
        {
            expression = "n = " + String.valueOf( n );
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "q = " + String.valueOf(q);
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "delta = 0.001;";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "m1 = Ceiling[ (1 + delta)*n*Log2[q] ]";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "r = " + String.valueOf(r);
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "l = Ceiling[Log[r,q]]";
            ml.evaluate(expression);
            ml.waitForAnswer();
            l = ml.getInteger();
            
            expression = "m2 = m1*l;";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "m = m1 + m2";
            ml.evaluate(expression);
            ml.waitForAnswer();
            this.m = ml.getInteger();
            SystemParameters.m = this.m;
            String out = "\n >> m : " + String.valueOf( this.m );
            
            expression = "A1nm = RandomChoice[Range[0,q-1],{n,m1}]";
            ml.evaluate(expression);
            ml.waitForAnswer();
            long[][] lA1nm = ml.getLongArray2();
            
            expression = "ns = NullSpace[A1nm, Modulus -> q];";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "nsI = Join[ns, q*IdentityMatrix[m1]];";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "{u, H} = HermiteDecomposition[nsI];";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "h2 = H[[1;;m1, 1;;m1]];";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "Transpose[h2]";
            ml.evaluate(expression);
            ml.waitForAnswer();
            long[][] lorto = ml.getLongArray2();
            
            LMatrix A1nm = new LMatrix(lA1nm);
            LMatrix orto = new LMatrix(lorto);
            
            AlwenPeikert ap = new AlwenPeikert( A1nm, orto, r, l, q );
            A = ap.setup( A1nm, folder, instanceName );
            
        }catch(MathLinkException mle)
        {
            System.err.println(" Error : MathLinkException computing TrapGen() " + 
                    mle.getMessage() );
        }
        
        return A;
    }

    private long Mod(long a, long q)
    {
        long mod;
        mod = a % q;
        mod = mod >= 0 ? mod : q + mod;
        
        return mod;
    }

    public LVector[] MultiExtract(PublicParam PP, LMatrix MK, LVector id, String folder)
    {
        System.out.println(" *********************  Extract  *****************");
        
        LMatrix HidB = H(id, 1);
        HidB = HidB.ModMult(PP.B, q);
        HidB = PP.A1.addition(HidB);
        
        LVector[] privateKeys = new LVector[PP.target.length];
        
        System.out.println(" computing trapdoor Gram-Shmidt ... \n");        
        String file = folder + "\\trapdoor.gso";
        
        MK.Compute_Gram_Shmidt();
        IO.saveGSO(file, MK.gso);
        //MK.printGSO();
        double cota = PP.sigma.doubleValue() * Math.sqrt( 2*m );
        System.out.println(" cota : " + cota);
        LVector e;
        for(int i=0; i < PP.target.length; i++)
        {
            do
            {
                e = Sample.SampleLeft2( PP, HidB, MK, PP.target[i], PP.e2[i]);
            }while( e.norm() > cota );
            System.out.println("|| e || : " + e.norm());
            privateKeys[i] = e;
        }
        
        return privateKeys;
    }
    
    
    public LVector[] MultiExtract_gso(PublicParam PP, LMatrix MK, LVector id, String gsopath)
    {
        System.out.println(" *********************  Extract  *****************");
        LMatrix HidB = H(id, 1);
        HidB = HidB.ModMult(PP.B, q);
        HidB = PP.A1.addition(HidB);
        
        LVector[] privateKeys = new LVector[PP.target.length];
        
        System.out.println(" reading trapdoor Gram-Shmidt ... \n");
        
        MK.gso = IO.readGSO( gsopath, MK.rows );
        
        System.out.println(" toBasis ...");
        MK.toBasis();

        double cota = PP.sigma.doubleValue() * Math.sqrt( 2*m );
        System.out.println(" cota : " + cota);
        LVector e;
        for(int i=0; i < PP.target.length; i++)
        {
            do
            {
                e = Sample.SampleLeft2( PP, HidB, MK, PP.target[i], PP.e2[i]);
            }while( e.norm() > cota );
            System.out.println("|| e || : " + e.norm());
            privateKeys[i] = e;
        }
        
        return privateKeys;
    }
    
    
    public CT Encrypt(PublicParam PP, LVector id, int b)
    {
        //System.out.println(" ***************** encrypt ****************");
        LMatrix aux;
        if( Fid == null)
        {            
            aux = H(id, 0).mult(PP.B);
            aux = PP.A1.addition(aux);
            
            Fid = PP.A0.augment(aux);
            
            Fid.reduceModule(q);
            
            FidT = Fid.transpose();
        }
        
        LVector s = RandomGen.Random_Zq_n();
        //System.out.print(" s : ");
        //s.print(); System.out.println("");
        
        LMatrix R = RandomGen.Random_R();
        /*System.out.println(" R : ");
        R.print(); System.out.println("\n\n");*/
        
        long x = RandomGen.Gaussian_Zq();
        //long x = 0; // noise zero
        //System.out.println(" x : " + x );
        
        LVector y = RandomGen.Gaussian_Zq_m();
        //LVector y = new LVector(m);  // noise zero
        //System.out.print(" y : ");
        //y.print(); System.out.println("");
        
        LMatrix Rt = R.transpose();
        LVector z = Rt.ModMult(y);
        //System.out.print(" z : ");
        //z.print(); System.out.println("");
        
        long c0 = PP.u[0].punto(s) + x + b*((long)Math.floor(q/2));
        c0 = Mod(c0, q);
        
        //aux = Fid.transpose();
        
        LVector vaux = FidT.ModMult(s);
        
        LVector c1 = vaux.add(y.augment(z));
        c1.reduceModule(q);
        
        return new CT(c0, c1);
    }
    
    public MultiCT MultiEncrypt(PublicParam PP, LVector id, byte[] block)
    {
        LMatrix aux;
        if( Fid == null)
        {            
            aux = H(id, 0).mult(PP.B);
            aux = PP.A1.addition(aux);
            
            Fid = PP.A0.augment(aux);
            
            Fid.reduceModule(q);
            
            FidT = Fid.transpose();
        }
        
        LVector s = RandomGen.Random_Zq_n();
        
        LMatrix R = RandomGen.Random_R();
        
        long x = 0; //RandomGen.Gaussian_Zq();
        
        long[] iy = new long[SystemParameters.m];
        LVector y = new LVector(iy); 
        //LVector y = RandomGen.Gaussian_Zq_m();
        
        LMatrix Rt = R.transpose();
        LVector z = Rt.ModMult(y);    
        
        byte[] block4bytes = new byte[4];
        System.arraycopy(block, 0, block4bytes, 1, 3);
        
        long ci;
        long[] c = new long[PP.u.length];
        int i, j, b, intblock = ByteBuffer.wrap(block4bytes).getInt();
        
        
        for( i=PP.u.length-1, j=0; i>=0; i--, j++)
        {
            b = (intblock >> i) & 1;
            ci = PP.u[j].punto(s) + x + b*((long)Math.floor(q/2));
            ci = Mod(ci, q);
            c[j] = ci;
            System.out.print(b);
        }
        
        LVector vaux = FidT.ModMult(s);
        
        LVector c1 = vaux.add( y.augment(z) );
        c1.reduceModule(q);
        
        return new MultiCT(c, c1);
    }
    
    public int Decrypt(PublicParam PP, LVector eid, CT ct)
    {
        long eidc1 = eid.punto(ct.c1);
        long w = ct.c0 - eidc1;
        long q2 = (long)Math.floor(q/2);
        long q4 = (long)Math.floor(q/4);
        long wq2;
        
        w = Mod(w, q);
        
        wq2 = Math.abs( w - q2 );
        
        /*
        System.out.println( " eidc1 : " + eidc1 );
        System.out.println( " c0    : " + ct.c0 );
        System.out.println( " w     : "  + w  );
        System.out.println( " q/2   : "  + q2 );
        System.out.println( " q/4   :  " + q4 );
        System.out.println( " w-q/2 : " + wq2 );
        */
        if( wq2 < q4 )
            return 1;
        else
            return 0;
    }

    public byte[] MultiDecrypt(PublicParam PP, PrivateKey prvkey, MultiCT ct)
     {
        LVector[] eid = prvkey.getLVectorArray();
        long eidc1;
        long q2 = (long)Math.floor(q/2);
        long q4 = (long)Math.floor(q/4);
        long w, wq2;
        
        ByteBuffer buff = ByteBuffer.allocate(4);
        
        int decrypted=0;
        
        for(int i=0; i < PP.u.length; i++)
        {
            eidc1 = eid[i].punto(ct.c1);
            w = ct.c0[i] - eidc1;
            w = Mod(w, q);
            
            wq2 = Math.abs( w - q2 );

            if( wq2 < q4 )
            {
                decrypted = (decrypted << 1) | 1;
                System.out.print("1");
            }
            else
            {
                decrypted = decrypted << 1;
                System.out.print("0");
            }
        }
        
        buff.putInt(decrypted);
        byte[] decryptedbyte = buff.array();
        byte[] b = new byte[3];
        System.arraycopy(decryptedbyte, 1, b, 0, 3);
        
        return b;
    }
    
    public LMatrix getSecret()
    {
        return Ta;
    }
    
    public static LMatrix H(LVector id, int encrypt_extract) // enc = 0 ; ext = 1
    {
        KernelLink ml = null;
        
        System.out.println(" Computing H(id) ....");
         //-------------------------------------------------------- 
         
        
            try 
            {
                ml = MathLinkFactory.createKernelLink("-linkmode launch -linkname "
                        + "'C:\\Program Files\\Wolfram Research\\Mathematica\\9.0\\MathKernel.exe'");
                if(SystemParameters.ml == null)
                    SystemParameters.ml = ml;
            }catch(MathLinkException mle)
            {
                System.out.println("Fatal error opening link : " + mle.getMessage());
                if(ml != null)
                    ml.close();
            }
        
        //--------------------------------------------------------        
        
        int i, j;
        
        String exp, Hstring, token1;
        long[][] H = new long[n][n];
        try
        {
            //if(encrypt_extract == 0)
            ml.discardAnswer();
            exp = "g[u_,x_] := Sum[u[[i]]*x^(i-1),{i,1,Length[u]}];";
            ml.evaluate(exp);
            ml.discardAnswer();
            
            exp = "u = " + id.getString();
            ml.evaluate(exp);
            ml.discardAnswer();
            
            exp = "f[i_] := CoefficientList[PolynomialMod[x^(i)*g[u,x]," + SystemParameters.irreduciblePoly + "],x]";
            ml.evaluate(exp);
            ml.discardAnswer();
            
            exp = "ToString[Array[f, " + String.valueOf(SystemParameters.n) + ",0]]";
            ml.evaluate(exp);
            ml.waitForAnswer();
            Hstring = ml.getString();
            
            // eliminate "{{"  and  "}}"
            Hstring = Hstring.substring(2,Hstring.length()-2);
            Hstring = Hstring.replaceAll("\\}, \\{", ":");
            StringTokenizer strTok1 = new StringTokenizer(Hstring, ":");
            StringTokenizer strTok2;
            
            i=0;
            while( strTok1.hasMoreTokens() )
            {
                token1 = strTok1.nextToken();
                strTok2 = new StringTokenizer(token1, ", ");
                
                j=0;
                while( strTok2.hasMoreTokens() )
                {
                    H[i][j] = Long.parseLong( strTok2.nextToken() );
                    j++;
                }
                i++;
            }
            
        }catch(MathLinkException mle)
        {
            System.err.println(" Error : MathLinkException computing H function " + 
                    mle.getMessage() );
        }
        
        return new LMatrix(H);
    } 
    
}