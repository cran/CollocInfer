SEIRRV = function(t,y,p,more)
{
    params = more$p.fun(t,more$p.more)
    beta = more$beta.fun(t,p,more$betadef)

    r = y
    names(r) = names(y)
    
    r[,'S'] = params[,'mu']*params[,'N']*(1-p['p']) + p['alpha3']*y[,'R2'] + p['phi']*params[,'gamma']*y[,'I'] - 
                y[,'S']*(p['i']+beta*y[,'I']) - params[,'mu']*y[,'S'] 
    r[,'E'] = y[,'S']*(p['i']+beta*y[,'I']) - params[,'sigma']*y[,'E'] - params[,'mu']*y[,'E']
    r[,'I'] = params[,'sigma']*y[,'E'] - params[,'gamma']*y[,'I'] - params[,'mu']*y[,'I']
    r[,'R1'] = (1-p['phi'])*params[,'gamma']*y[,'I'] + p['xi']*(p['i']+beta*y[,'I'])*y[,'R2'] - params[,'mu']*y[,'R1'] - p['alpha1']*y[,'R1'] 
    r[,'R2'] = p['alpha1']*y[,'R1'] + p['alpha2']*y[,'V'] - p['xi']*(p['i']+beta*y[,'I'])*y[,'R2'] - p['alpha3']*y[,'R1']
    r[,'V'] = params[,'mu']*(1-p['p'])*params[,'N'] - p['alpha2']*y[,'V'] - params[,'mu']*y[,'V']

    return(r)
}


SEIRV = function(t,y,p,more)
{
    params = more$p.fun(t,more$p.more)
    beta = more$beta.fun(t,p,more$betadef)

    r = y
    names(r) = names(y)
    
    r[,'S'] = params[,'mu']*params[,'N']*(1-p['p']) + p['alpha']*y[,'R'] + p['phi']*params[,'gamma']*y[,'I'] + 
                p['alpha2']*y[,'V'] - y[,'S']*(p['i']+beta*y[,'I']) - params[,'mu']*y[,'S'] 
    r[,'E'] = y[,'S']*(p['i']+beta*y[,'I']) - params[,'sigma']*y[,'E'] - params[,'mu']*y[,'E']
    r[,'I'] = params[,'sigma']*y[,'E'] - params[,'gamma']*y[,'I'] - params[,'mu']*y[,'I']
    r[,'R'] = (1-p['phi'])*params[,'gamma']*y[,'I']  - params[,'mu']*y[,'R'] - p['alpha']*y[,'R'] 
    r[,'V'] = params[,'mu']*(1-p['p'])*params[,'N'] - p['alpha2']*y[,'V'] - params[,'mu']*y[,'V']

    return(r)
}



SEIRS = function(t,y,p,more)
{
    params = more$p.fun(t,more$p.more)
    beta = more$beta.fun(t,p,more$betadef)

    r = y
    names(r) = names(y)
    
    r[,'S'] = params[,'mu']*params[,'N']*(1-p['p']) + p['alpha']*y[,'R'] + 
                 - y[,'S']*(p['i']+beta*y[,'I']) - params[,'mu']*y[,'S'] 
    r[,'E'] = y[,'S']*(p['i']+beta*y[,'I']) - params[,'sigma']*y[,'E'] - params[,'mu']*y[,'E']
    r[,'I'] = params[,'sigma']*y[,'E'] - params[,'gamma']*y[,'I'] - params[,'mu']*y[,'I']
    r[,'R'] = params[,'gamma']*y[,'I'] - params[,'mu']*y[,'R'] - p['alpha']*y[,'R'] 

    return(r)
}


SEIR = function(t,y,p,more)
{
    params = more$p.fun(t,more$p.more)
    beta = more$beta.fun(t,p,more$betadef)

    r = y
    names(r) = names(y)

    r[,'S'] = params[,'mu']*params[,'N']*(1-p['p']) +
                 - y[,'S']*(p['i']+beta*y[,'I']) - params[,'mu']*y[,'S']
    r[,'E'] = y[,'S']*(p['i']+beta*y[,'I']) - params[,'sigma']*y[,'E'] - params[,'mu']*y[,'E']
    r[,'I'] = params[,'sigma']*y[,'E'] - params[,'gamma']*y[,'I'] - params[,'mu']*y[,'I']
    r[,'R'] = params[,'gamma']*y[,'I'] - params[,'mu']*y[,'R']

    return(r)
}