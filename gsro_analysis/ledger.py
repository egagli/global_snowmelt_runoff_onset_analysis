"""The icechunk ledger shared by the analysis stores (the ERA5-Land store and
the ancillary store): Azure storage handles, the repository config with
storage retries for many concurrent GitHub Actions writers, the bounded
commit-with-retry that rebases over concurrent commits, and the fold over the
commit history that turns metadata-carrying commits into a work ledger — the
production repo's fleet pattern (docs/icechunk-github-actions-pattern.md there).

A work unit (a water year, a tile) writes disjoint chunks and makes ONE commit
with machine-readable metadata; a failed unit commits nothing; "what remains"
is what has no commit. Commit messages are for humans, tooling reads only the
metadata (``schema`` versions it).
"""

import random
import time

SCHEMA = 1
BRANCH = "main"
COMMIT_MAX_TRIES = 6


def azure_storage(config, prefix):
    """icechunk storage for a container-qualified Azure prefix ('snowmelt/...')."""
    import icechunk
    container, path = prefix.split('/', 1)
    return icechunk.azure_storage(account=config.azure_storage_account, container=container,
                                  prefix=path, sas_token=config.sas_token)


def storage(config, prefix, local_store=None):
    """``local_store`` (a directory path) selects a local filesystem repo for tests."""
    import icechunk
    return icechunk.local_filesystem_storage(str(local_store)) if local_store else azure_storage(config, prefix)


def repo_config():
    """Storage retries for many concurrent writers (as the production output repo)."""
    import icechunk
    cfg = icechunk.RepositoryConfig.default()
    cfg.storage = icechunk.StorageSettings(retries=icechunk.StorageRetriesSettings(
        max_tries=20, initial_backoff_ms=200, max_backoff_ms=60_000))
    return cfg


def repo_exists(config, prefix, local_store=None):
    import icechunk
    return icechunk.Repository.exists(storage(config, prefix, local_store))


def open_repo(config, prefix, local_store=None, what='repository'):
    """Open a repository; a clear FileNotFoundError when it has not been created yet."""
    import icechunk
    st = storage(config, prefix, local_store)
    if not icechunk.Repository.exists(st):
        raise FileNotFoundError(f"no {what} at {local_store or prefix}: its workflow creates it")
    return icechunk.Repository.open(st, config=repo_config())


def create_repo(config, prefix, local_store=None):
    import icechunk
    repo = icechunk.Repository.create(storage(config, prefix, local_store), config=repo_config())
    repo.save_config()
    return repo


def delete_repo(config, prefix, local_store=None, log=print):
    """The only deletion in this module — the workflows' off-by-default start_fresh boxes."""
    if local_store:
        import shutil
        shutil.rmtree(local_store, ignore_errors=True)
    else:
        from gsro_analysis import settings
        fs = settings.fresh_blob_fs(config)
        if fs.exists(prefix):
            fs.rm(prefix, recursive=True)
        fs.invalidate_cache()
    log(f"deleted {local_store or prefix} (start_fresh)")


def provenance():
    try:
        from global_snowmelt_runoff_onset.provenance import collect_provenance
        return collect_provenance()
    except Exception:  # noqa: BLE001 - provenance must never fail a commit
        import platform
        return {"host": platform.node()}


def commit_records(repo, branch=BRANCH):
    """Newest -> oldest list of the pipeline commits' metadata (commits without the
    schema key — init, maintenance — are skipped), each with ``ancestry_index``
    (0 = newest), ``snapshot_id`` and ``written_at``."""
    records = []
    for index, snap in enumerate(repo.ancestry(branch=branch)):
        meta = snap.metadata or {}
        if 'schema' not in meta or 'kind' not in meta:
            continue
        records.append({'ancestry_index': index, 'snapshot_id': snap.id,
                        'written_at': str(snap.written_at), **meta})
    return records


def commit_with_retry(repo, write_fn, message, metadata, branch=BRANCH, allow_empty=False,
                      max_tries=COMMIT_MAX_TRIES, log=print):
    """Write through a FRESH session and commit, rebasing over concurrent commits
    (``ConflictDetector``: always clean here, the writers touch disjoint chunks) and
    retrying the whole write on transient errors — the production fleet's routine.
    Programming/schema errors re-raise immediately."""
    import icechunk
    last = None
    for attempt in range(max_tries):
        try:
            session = repo.writable_session(branch)
            write_fn(session)
            return session.commit(message, metadata=metadata,
                                  rebase_with=icechunk.ConflictDetector(), allow_empty=allow_empty)
        except (ValueError, KeyError, TypeError, AssertionError):
            raise
        except Exception as e:  # noqa: BLE001 - conflict, expired session, storage blip
            last = e
            delay = min(60, 2 ** attempt) * random.uniform(0.5, 1.5)
            log(f"commit attempt {attempt + 1}/{max_tries} failed ({type(e).__name__}: {e}); retry in {delay:.0f}s")
            time.sleep(delay)
    raise RuntimeError(f"commit failed after {max_tries} attempts") from last
