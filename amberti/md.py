from amberti.amber import pmemd
from amberti.logger import getLogger
import os

logger = getLogger()


def em(
        defname,
        prmtop, 
        conf,
        ref,
        maxcyc=10000,
        resstraint_wt=5.00,
        fep=True,
        clambda=None,
        scalpha=None,
        scbeta=None,
        ifsc=1,
        crgmask=None,
        timask1=None,
        timask2=None,
        scmask1=None,
        scmask2=None,
        fname="min.in",
        run=True,
    ):

    script = [
        defname,
        "&cntrl",
        "imin = 1, ntmin = 2,",
        f"maxcyc = {maxcyc},",
        "ntpr = 20, ntwe = 20,",
        "ntb = 1,",
        f"ntr = 1, restraint_wt = {resstraint_wt},",
        "restraintmask='!:WAT & !@H=',",
    ]
    if fep:
        script += [
            f"icfe = 1, clambda = {clambda}, scalpha = {scalpha}, scbeta = {scbeta},",
            "logdvdl = 0,",
            f"timask1 = '{timask1}', timask2 = '{timask2}',",
            f"ifsc = {ifsc},"
        ]
        if crgmask is not None:
            script.append(f"crgmask={crgmask},")
        if ifsc == 1:
            assert scmask1 is not None, "you are setting ifsc=1 indicates that soft core" \
            "is activated, however, scmask1 is not assigned."
            assert scmask2 is not None, "you are setting ifsc=1 indicates that soft core" \
            "is activated, however, scmask2 is not assigned."
            script.append(f"scmask1 = '{scmask1}', scmask2 = '{scmask2}',")

    script += [
        "/",
        "&ewald",
        "/"
    ]

    with open(fname, "w") as fp:
        fp.write("\n".join(script))
        fp.write("\n")
    
    return_code = pmemd(defname, fname, prmtop, conf, outtraj=False, ref=ref, cuda=True, mpi=False, run=run)
    return return_code


def heat(
        defname,
        prmtop, 
        conf,
        ref,
        nsteps=10000,
        dt=0.002,
        temp=300,
        resstraint_wt=5.00,
        irest=0, ntx=1,
        fep=True,
        clambda=None,
        scalpha=None,
        scbeta=None,
        ifsc=1,
        crgmask=None,
        timask1=None,
        timask2=None,
        scmask1=None,
        scmask2=None,
        tempi=5.0,
        ofreq=1000,
        fname="heat.in",
        run=True,
    ):
    
    script = [
        defname,
        "&cntrl",
        f"imin = 0, nstlim = {nsteps}, irest = {irest}, ntx = {ntx}, dt = {dt},",
        "nmropt = 1,",
        f"ntt = 1, temp0 = {temp}, tempi = {tempi}, tautp = 1.0,",
        "ntb = 1,",
        "ntc = 2, ntf = 1,",
        "ioutfm = 1, iwrap = 1,",
        f"ntwe = {ofreq}, ntwx = {ofreq}, ntpr = {ofreq}, ntwr = {ofreq},",
        f"ntr = 1, restraint_wt = {resstraint_wt},",
        "restraintmask='!:WAT & !@H=',",
    ]

    if fep:
        script += [
            f"icfe = 1, clambda = {clambda}, scalpha = {scalpha}, scbeta = {scbeta},",
            "logdvdl = 0,",
            f"timask1 = '{timask1}', timask2 = '{timask2}',",
            f"ifsc = {ifsc},"
        ]
        if crgmask is not None:
            script.append(f"crgmask={crgmask},")
        if ifsc == 1:
            assert scmask1 is not None, "you are setting ifsc=1 indicates that soft core" \
            "is activated, however, scmask1 is not assigned."
            assert scmask2 is not None, "you are setting ifsc=1 indicates that soft core" \
            "is activated, however, scmask2 is not assigned."
            script.append(f"scmask1 = '{scmask1}', scmask2 = '{scmask2}',")
    
    script += [
        "/",
        "&ewald",
        "/" ,
        "&wt",
        "type='TEMP0',",
        f"istep1 = 0, istep2 = {int(nsteps*0.8)},",                         
        f"value1 = {tempi}, value2 = {temp}",
        "/",
        "&wt type = 'END'",
        "/",
    ]
    
    with open(fname, "w") as fp:
        fp.write("\n".join(script))
        fp.write("\n")
    
    return_code = pmemd(defname, fname, prmtop, conf, outtraj=True, ref=ref, cuda=True, mpi=False, run=run)
    return return_code


def pressurize(
        defname,
        prmtop, 
        conf,
        ref,
        nsteps=10000,
        dt=0.002,
        temp=300,
        resstraint_wt=5.00,
        irest=0, ntx=1,
        fep=True,
        clambda=None,
        scalpha=None,
        scbeta=None,
        ifsc=1,
        crgmask=None,
        timask1=None,
        timask2=None,
        scmask1=None,
        scmask2=None,
        ofreq=1000,
        fname="press.in",
        run=True
    ):

    script = [
        defname,
        "&cntrl",
        f"imin = 0, nstlim = {nsteps}, irest = {irest}, ntx = {ntx}, dt = {dt},",
        f"ntt = 1, temp0 = {temp}, tautp = 1.0,",
        "ntp = 1, pres0 = 1.0, taup = 2.0,",
        "ntb = 2,",
        "ntc = 2, ntf = 1,",
        "ioutfm = 1, iwrap = 1,",
        f"ntwe = {ofreq}, ntwx = {ofreq}, ntpr = {ofreq}, ntwr = {ofreq},",

        f"ntr = 1, restraint_wt = {resstraint_wt},",
        "restraintmask='!:WAT & !@H=',",
    ]

    if fep:
        script += [
            f"icfe = 1, clambda = {clambda}, scalpha = {scalpha}, scbeta = {scbeta},",
            "logdvdl = 0,",
            f"timask1 = '{timask1}', timask2 = '{timask2}',",
            f"ifsc = {ifsc},"
        ]
        if crgmask is not None:
            script.append(f"crgmask={crgmask},")
        if ifsc == 1:
            assert scmask1 is not None, "you are setting ifsc=1 indicates that soft core" \
            "is activated, however, scmask1 is not assigned."
            assert scmask2 is not None, "you are setting ifsc=1 indicates that soft core" \
            "is activated, however, scmask2 is not assigned."
            script.append(f"scmask1 = '{scmask1}', scmask2 = '{scmask2}',")

    script += [
        "/",
        "&ewald",
        "/" ,
    ]

    with open(fname, "w") as fp:
        fp.write("\n".join(script))
        fp.write("\n")
    
    return_code = pmemd(defname, fname, prmtop, conf, outtraj=True, ref=ref, cuda=True, mpi=False, run=run)
    return return_code


def equilibrium(
        prmtop, 
        conf,
        ref,
        maxcyc=10000,
        nsteps_heat=10000,
        nsteps_press=50000,
        temp=300,
        resstraint_wt=5.00,
        fep=True,
        timask1=None,
        timask2=None,
        scmask1=None,
        scmask2=None,
        ofreq=1000,
    ):

    logger.info("Emergy minimizing ...")
    em(
        "min",
        prmtop, 
        conf,
        ref,
        maxcyc=maxcyc,
        resstraint_wt=resstraint_wt,
        fep=True,
        timask1=timask1,
        timask2=timask2,
        scmask1=scmask1,
        scmask2=scmask2,
        fname="min.in"
    )
    assert os.path.exists("min.rst7")

    logger.info("Pre-heating ...")
    heat(
        defname="pre_heating",
        prmtop=prmtop, 
        conf="min.rst7",
        ref="min.rst7",
        nsteps=20,
        dt=0.0002,
        temp=5.0,
        resstraint_wt=resstraint_wt,
        fep=fep,
        timask1=timask1,
        timask2=timask2,
        scmask1=scmask1,
        scmask2=scmask2,
        tempi=4.0,
        ofreq=1,
        fname="pre_heat.in"
    )
    assert os.path.exists("pre_heat.rst7")

    logger.info("Second Emergy minimizing ...")
    em(
        "min2",
        prmtop=prmtop,
        conf="pre_heating.rst7",
        ref="pre_heating.rst7",
        maxcyc=maxcyc,
        resstraint_wt=resstraint_wt,
        fep=True,
        timask1=timask1,
        timask2=timask2,
        scmask1=scmask1,
        scmask2=scmask2,
        fname="min2.in"
    )
    assert os.path.exists("min2.rst7")

    logger.info(f"Heating from 5K to {temp}K ...")
    heat(
        defname="heating",
        prmtop=prmtop, 
        conf="min2.rst7",
        ref="min2.rst7",
        nsteps=nsteps_heat,
        temp=temp,
        resstraint_wt=resstraint_wt,
        fep=fep,
        timask1=timask1,
        timask2=timask2,
        scmask1=scmask1,
        scmask2=scmask2,
        tempi=5.0,
        ofreq=ofreq,
        fname="heat.in"
    )
    assert os.path.exists("heat.rst7")

    logger.info("Pre-Pressurising ...")
    pressurize(
        defname="pre_pressing",
        prmtop=prmtop, 
        conf="heat.rst7",
        ref="heat.rst7",
        nsteps=1000,
        dt=0.002,
        temp=temp,
        resstraint_wt=resstraint_wt,
        irest=1, ntx=5,
        fep=True,
        timask1=timask1,
        timask2=timask2,
        scmask1=scmask1,
        scmask2=scmask2,
        ofreq=10,
        fname="pre_press.in"
    )
    assert os.path.exists("pre_pressing.rst7")

    logger.info("Pressurising ...")
    pressurize(
        defname="pressing",
        prmtop=prmtop, 
        conf="pre_press.rst7",
        ref="pre_press.rst7",
        nsteps=nsteps_press,
        dt=0.002,
        temp=temp,
        resstraint_wt=resstraint_wt,
        irest=1, ntx=5,
        fep=True,
        timask1=timask1,
        timask2=timask2,
        scmask1=scmask1,
        scmask2=scmask2,
        ofreq=1000,
        fname="pre_press.in"
    )
    assert os.path.exists("pressing.rst7")

    logger.info("Equilibrium finished.")

    return 0
    

def production(
        defname,
        prmtop, 
        conf,
        ref,
        nsteps=10000,
        dt=0.002,
        temp=300,
        resstraint_wt=5.00,
        irest=1, ntx=5,
        fep=True,
        clambda=None,
        scalpha=None,
        scbeta=None,
        ifmbar=1,
        mbar_states=None,
        mbar_lambda=None,
        ifsc=0,
        crgmask=None,
        timask1=None,
        timask2=None,
        scmask1=None,
        scmask2=None,
        ntwe=1000,
        ntwx=10000,
        ntpr=10000,
        ntwr=20000,
        fname="press.in",
        run=False
    ):
    """Run production MD.
        ntwx:    Every ntwx steps, the coordinates will be written to the mdcrd file.
        ntwr:    Every ntwr steps during dynamics, the “restrt” file will be written.
        ntpr:    Every ntpr steps, energy information will be printed in human-readable form to
                 files "mdout" and "mdinfo".
        ntwe:    Every ntwe steps, the energies and temperatures will be written to file "mden" in a compact form.
        ioutfm:  The format of coordinate and velocity trajectory files. 
                    *  1 for binary NetCDF. 
                    *  0 for ASCII.
    """
    script = [
        defname,
        "&cntrl",
        f"imin = 0, nstlim = {nsteps}, irest = {irest}, ntx = {ntx}, dt = {dt},",
        f"ntt = 3, temp0 = {temp}, gamma_ln = 2.0, ig = -1,",
        "ntp = 1, pres0 = 1.0, taup = 2.0,",
        "ntb = 2,",
        "ntc = 2, ntf = 1,",
        "ioutfm = 1, iwrap = 1,",
        f"ntwe = {ntwe}, ntwx = {ntwx}, ntpr = {ntpr}, ntwr = {ntwr},",

        f"ntr = 1, restraint_wt = {resstraint_wt},",
        "restraintmask='!:WAT & !@H=',",
    ]

    if fep:
        script += [
            f"icfe = 1, clambda = {clambda}, scalpha = {scalpha}, scbeta = {scbeta},",
            "logdvdl = 0,",
            f"timask1 = '{timask1}', timask2 = '{timask2}',",
            f"ifsc = {ifsc}, crgmask = '{crgmask}',"
            
        ]
        if scmask1 is not None:
            script += [
                f"scmask1='{scmask1}', scmask2='{scmask2}',"
            ]
    
    if fep:
        script += [
            f"icfe = 1, clambda = {clambda}, scalpha = {scalpha}, scbeta = {scbeta},",
            "logdvdl = 0,",
            f"timask1 = '{timask1}', timask2 = '{timask2}',",
            f"ifsc = {ifsc}, ifmbar = {ifmbar},"
        ]
        if crgmask is not None:
            script.append(f"crgmask={crgmask},")
        if ifsc == 1:
            assert scmask1 is not None, "you are setting ifsc=1 indicates that soft core" \
            "is activated, however, scmask1 is not assigned."
            assert scmask2 is not None, "you are setting ifsc=1 indicates that soft core" \
            "is activated, however, scmask2 is not assigned."
            script.append(f"scmask1 = '{scmask1}', scmask2 = '{scmask2}',")
        if ifmbar == 1:
            assert mbar_states is not None, "you are setting ifmbar=1 indicates that MBAR" \
            "is activated, however, mbar_states is not assigned."
            assert mbar_lambda is not None, "you are setting ifmbar=1 indicates that MBAR" \
            "is activated, however, mbar_lambda is not assigned."
            script.append(f"mbar_states = {mbar_states}, mbar_lambda = {mbar_lambda},")

    script += [
        "/",
        "&ewald",
        "/" ,
    ]
    
    with open(fname, "w") as fp:
        fp.write("\n".join(script))
        fp.write("\n")

    return_code = pmemd(defname, fname, prmtop, conf, outtraj=True, ref=ref, cuda=True, mpi=False, run=run)
    return return_code


