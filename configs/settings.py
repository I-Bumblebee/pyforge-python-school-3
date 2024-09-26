from os import getenv

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    database_url: str
    redis_port: int
    redis_host: str
    cache_duration: int
    aws_access_key_id: str
    aws_secret_access_key: str
    region_name: str

    model_config = SettingsConfigDict(
        case_sensitive=False,
        env_file=(
            ".env.testing" if getenv("APP_ENV") == "testing" else ".env"
        ),
        env_file_encoding="utf-8",
        extra="ignore",
    )


settings = Settings()
