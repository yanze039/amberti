import logging
import os
import sys


logging.basicConfig(
    format=">> | %(asctime)s | %(levelname)s | %(name)s | %(message)s |",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)

def getLogger():
    # Setting logging configurations
    # logging.root.setLevel(0)
    logger = logging.getLogger(__name__)
    return logger