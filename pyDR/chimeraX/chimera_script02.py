import sys
from chimerax.core.commands import run as rc

rc(session,"remotecontrol rest start port 60960")
sys.path.append("/Users/albertsmith/Documents/GitHub/GHSR_archive/pyDR/chimeraX")
sys.path.append("/Users/albertsmith/Documents/GitHub/GHSR_archive/pyDR")
from RemoteCMXside import CMXReceiver as CMXR
import RemoteCMXside
cmxr=CMXR(session,7002,rc_port0=60960)
rc(session,"ui mousemode right select")
