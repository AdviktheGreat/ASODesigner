"""Tests for the logging utilities."""

import logging

from aso_designer.utils.logger import get_logger, setup_logger


class TestGetLogger:
    def test_returns_logger_with_handler(self):
        logger = get_logger("aso_designer.test_a")
        assert isinstance(logger, logging.Logger)
        assert logger.handlers  # console handler configured

    def test_repeated_calls_reuse_handlers(self):
        """A second call returns the same configured logger, no duplicates."""
        logger = get_logger("aso_designer.test_b")
        n_handlers = len(logger.handlers)
        again = get_logger("aso_designer.test_b")
        assert again is logger
        assert len(again.handlers) == n_handlers

    def test_default_name(self):
        logger = get_logger()
        assert logger.name == "aso_designer"


class TestSetupLogger:
    def test_sets_level(self):
        logger = setup_logger("aso_designer.test_c", level=logging.DEBUG)
        assert logger.level == logging.DEBUG

    def test_log_file_created(self, tmp_path):
        log_file = tmp_path / "logs" / "aso.log"
        logger = setup_logger("aso_designer.test_d", log_file=str(log_file))
        logger.info("hello file")
        for handler in logger.handlers:
            handler.flush()
        assert log_file.exists()

    def test_duplicate_handlers_not_added(self, tmp_path):
        logger = setup_logger("aso_designer.test_e")
        n = len(logger.handlers)
        setup_logger("aso_designer.test_e")
        assert len(logger.handlers) == n
