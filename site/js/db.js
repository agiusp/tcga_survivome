// DuckDB-WASM: loads the patient table once and answers cohort queries with SQL.
import * as duckdb from 'https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@1.29.0/+esm';

let conn;

export async function initDb(parquetUrl) {
  const bundle = await duckdb.selectBundle(duckdb.getJsDelivrBundles());
  // Workers can't be created from another origin directly, so wrap the CDN script.
  const workerUrl = URL.createObjectURL(
    new Blob([`importScripts("${bundle.mainWorker}");`], { type: 'text/javascript' }));
  const db = new duckdb.AsyncDuckDB(new duckdb.ConsoleLogger(duckdb.LogLevel.WARNING), new Worker(workerUrl));
  await db.instantiate(bundle.mainModule, bundle.pthreadWorker);
  URL.revokeObjectURL(workerUrl);
  const res = await fetch(parquetUrl);
  if (!res.ok) throw new Error(`Could not load ${parquetUrl} (${res.status})`);
  await db.registerFileBuffer('survivome.parquet', new Uint8Array(await res.arrayBuffer()));
  conn = await db.connect();
  await conn.query(`CREATE TABLE s AS SELECT * FROM parquet_scan('survivome.parquet')`);
}

/** Run SQL; returns plain row objects (BigInt counts become numbers). */
export async function query(sql) {
  const result = await conn.query(sql);
  return result.toArray().map((row) => {
    const o = row.toJSON();
    for (const k in o) if (typeof o[k] === 'bigint') o[k] = Number(o[k]);
    return o;
  });
}

export const ident = (name) => `"${String(name).replace(/"/g, '""')}"`;
export const lit = (value) => `'${String(value).replace(/'/g, "''")}'`;
