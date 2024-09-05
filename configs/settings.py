from os import getenv

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    database_url: str
    redis_port: int
    redis_host: str
    cache_duration: int

    model_config = SettingsConfigDict(
        case_sensitive=False,
        env_file=(
            ".env.testing" if getenv("APP_ENV") == "testing" else ".env"
        ),
        env_file_encoding="utf-8",
        extra="ignore",
    )


settings = Settings()
