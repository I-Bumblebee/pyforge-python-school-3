import logging.config
from typing import AsyncGenerator

import pytest_asyncio
from httpx import ASGITransport, AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession, create_async_engine
from sqlalchemy.orm import sessionmaker

from configs.database import get_db_session
from configs.settings import settings
from configs.setup_logger import setup_logger
from src.main import app

logger = logging.getLogger(__name__)
setup_logger()


@pytest_asyncio.fixture(scope="session")
async def async_client() -> AsyncGenerator[AsyncClient, None]:
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as ac:
        yield ac


@pytest_asyncio.fixture(scope="function")
async def test_session() -> AsyncGenerator[AsyncSession, None]:
    AsyncSessionTesting = sessionmaker(
        bind=create_async_engine(settings.database_url),
        class_=AsyncSession,
        autoflush=True,
    )
    async with AsyncSessionTesting() as session:
        async with session.begin():
            app.dependency_overrides[get_db_session] = lambda: session

            yield session

            await session.rollback()

            app.dependency_overrides.pop(get_db_session, None)
