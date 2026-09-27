"""Structured logging configuration for ASODesigner.

Provides a configurable logger with console and rotating file handlers.
Inspired by GuideX's basic logger but extended with file output, log
rotation, and configurable verbosity levels.
"""

import logging
import sys
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Optional, Union

DEFAULT_FORMAT = "%(asctime)s | %(name)s | %(levelname)s | %(message)s"
DEFAULT_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def setup_logger(
    name: str = "asodesigner",
    level: Union[int, str] = logging.INFO,
    log_file: Optional[Union[str, Path]] = None,
    max_bytes: int = 5 * 1024 * 1024,
    backup_count: int = 3,
) -> logging.Logger:
    """Configure and return a logger with console and optional file handlers.

    Parameters
    ----------
    name : str
        Logger name (typically module name).
    level : int or str
        Logging level (e.g. logging.DEBUG, "INFO").
    log_file : str or Path, optional
        Path to a log file. If provided, a rotating file handler is attached.
    max_bytes : int
        Maximum file size before rotation (default 5 MB).
    backup_count : int
        Number of rotated backup files to keep.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger(name)

    # Avoid duplicate handlers on repeated calls
    if logger.handlers:
        return logger

    logger.setLevel(level)

    formatter = logging.Formatter(DEFAULT_FORMAT, datefmt=DEFAULT_DATE_FORMAT)

    # Console handler — always present
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Rotating file handler — optional
    if log_file is not None:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = RotatingFileHandler(
            log_path,
            maxBytes=max_bytes,
            backupCount=backup_count,
            encoding="utf-8",
        )
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    logger.propagate = False
    return logger


def get_logger(name: str = "asodesigner") -> logging.Logger:
    """Return an existing logger or create a default one."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        setup_logger(name)
    return logger
