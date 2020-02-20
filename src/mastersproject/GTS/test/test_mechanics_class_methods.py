import logging
from typing import (  # noqa
    Any,
    Coroutine,
    Generator,
    Generic,
    Iterable,
    List,
    Mapping,
    Optional,
    Set,
    Tuple,
    Type,
    TypeVar,
    Union,
)
import os
from pathlib import Path

import porepy as pp
import numpy as np
from porepy.models.contact_mechanics_model import ContactMechanics

import GTS as gts
from util.logging_util import __setup_logging
logger = logging.getLogger(__name__)

