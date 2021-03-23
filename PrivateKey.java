package ABB;

import java.io.Serializable;

public class PrivateKey implements Serializable
{
    private LVector[] prvkey;
    
    public PrivateKey(LVector[] prvkey)
    {
        this.prvkey = prvkey;
    }
    
    public LVector[] getLVectorArray()
    {
        return this.prvkey;
    }
}
