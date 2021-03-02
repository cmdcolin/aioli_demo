import { useCallback, useEffect, useRef, useState } from "react";
import { Aioli } from "@biowasm/aioli";
import GranularRectLayout from "./layout";
import { StringParam, useQueryParams, withDefault } from "use-query-params";

function parseLocString(locString) {
  const [refId, rest] = locString.split(":");
  const [start, end] = rest.split("-");
  return { refId: +refId - 1, start: +start, end: +end };
}

export function parseCigar(cigar) {
  return (cigar || "").split(/([MIDNSHPX=])/);
}

// because we are using use-query-params without a router
export function useForceUpdate() {
  const [, setTick] = useState(0);
  const update = useCallback(() => {
    setTick((tick) => tick + 1);
  }, []);
  return update;
}

const featHeight = 10;
const width = 1800;
const height = 1000;
const snpcovheight = 100;
const initialLoc = "1:20000-40000";
const initialFile =
  "https://s3.amazonaws.com/1000genomes/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam";
const initialFasta = "https://jbrowse.org/genomes/hg19/fasta/hg19.fa.gz";

function App() {
  const ref = useRef();
  const snpcovref = useRef();
  const [params, setParams] = useQueryParams({
    loc: withDefault(StringParam, initialLoc),
    file: withDefault(StringParam, initialFile),
    fasta: withDefault(StringParam, initialFasta),
  });
  const [readData, setReadData] = useState();
  const [mpileupData, setMPileupData] = useState();
  const [samtools, setSamtools] = useState();
  const [bamFile, setBamFile] = useState();
  const [fastaFile, setFastaFile] = useState();
  const [file, setFile] = useState(params.file);
  const [loc, setLoc] = useState(params.loc);
  const forceUpdate = useForceUpdate();

  // this block initializes the samtools tool on Aioli
  useEffect(() => {
    (async () => {
      const tool = new Aioli("samtools/1.10");
      await tool.init();
      setSamtools(tool);
    })();
  }, []);

  // this block "mounts" the BAM file on the Aioli FS
  useEffect(() => {
    (async () => {
      // wait on "samtools" being initialized even though "samtools" is not
      // used here, because it must mount something internally in the
      // Aioli/samtools interface
      if (samtools) {
        const url = new URL(params.file, window.location);
        const fasta = await Aioli.mount(`${params.fasta}`);
        const bam = await Aioli.mount(`${url}`);
        await Aioli.mount(`${url}.bai`);

        console.log({ fasta });
        setBamFile(bam);
        setFastaFile(fasta);
      }
    })();
  }, [params.file, samtools, params.fasta]);

  // this block performs a `samtools view` query
  useEffect(() => {
    (async () => {
      if (bamFile && samtools) {
        const d = await samtools.exec(`view ${bamFile.path} ${params.loc}`);
        setReadData(d);
      }
    })();
  }, [bamFile, params.loc, samtools]);

  // this block performs a `samtools mpileup` query
  useEffect(() => {
    (async () => {
      if (bamFile && samtools && fastaFile) {
        const d = await samtools.exec(
          `mpileup -r ${params.loc} ${bamFile.path}`
        );
        setMPileupData(d);
      }
    })();
  }, [bamFile, params.loc, samtools, fastaFile]);

  // this block draws the rectangles
  useEffect(() => {
    const ctx = ref.current.getContext("2d");
    ctx.clearRect(0, 0, width, height);
    const parsedLoc = parseLocString(params.loc);
    const bpPerPx = width / (parsedLoc.end - parsedLoc.start);
    const layout = new GranularRectLayout();
    readData?.stdout.split("\n").forEach((row, index) => {
      const [, flagString, , startString, , CIGAR] = row.split("\t");
      const start = +startString;
      const flags = +flagString;
      const cigarOps = parseCigar(CIGAR);
      let length = 0;
      for (let i = 0; i < cigarOps.length; i += 2) {
        const len = +cigarOps[i];
        const op = cigarOps[i + 1];
        if (op !== "I" && op !== "S" && op !== "H") {
          length += len;
        }
      }
      if (flags & 16) {
        ctx.fillStyle = "#f99";
      } else {
        ctx.fillStyle = "#99f";
      }

      const end = start + length;
      const leftPx = (start - parsedLoc.start) * bpPerPx;
      const width = (end - start) * bpPerPx;
      const y = layout.addRect(index, start, end, featHeight);

      ctx.fillRect(leftPx, y, width, featHeight);
    });
  }, [readData, params.loc]);

  // this block draws the rectangles
  useEffect(() => {
    const ctx = snpcovref.current.getContext("2d");
    ctx.clearRect(0, 0, width, snpcovheight);
    const parsedLoc = parseLocString(params.loc);
    const bpPerPx = width / (parsedLoc.end - parsedLoc.start);
    mpileupData?.stdout.split("\n").forEach((row) => {
      const [, startString, , numReadsString, bases] = row.split("\t");
      const start = +startString;
      const numReads = +numReadsString;
      const leftPx = (start - parsedLoc.start) * bpPerPx;
      const end = start + 1;
      const width = (end - start) * bpPerPx;

      ctx.fillStyle = "#ccc";
      ctx.fillRect(leftPx, snpcovheight - numReads, width, numReads);
      const map = {};
      bases?.split("").forEach((base) => {
        const b = base.toUpperCase();
        map[b] = (map[b] || 0) + 1;
      });
      let curr = 0;
      const colors = { A: "green", C: "blue", G: "yellow", T: "red" };
      console.log({ map });
      ["A", "C", "G", "T"].forEach((base) => {
        curr += map[base];
        ctx.fillStyle = colors[base];
        ctx.fillRect(leftPx, snpcovheight - curr, 1, map[base]);
      });
    });
  }, [mpileupData, params.loc]);

  return (
    <div>
      <p>Enter a BAM file URL:</p>
      <form
        onSubmit={(event) => {
          setParams({ file, loc });
          setBamFile(undefined);
          forceUpdate();
          event.preventDefault();
        }}
      >
        <label htmlFor="url" />
        <input
          id="url"
          type="text"
          value={file}
          onChange={(event) => setFile(event.target.value)}
        />

        <label htmlFor="loc" />
        <input
          id="loc"
          type="text"
          value={loc}
          onChange={(event) => setLoc(event.target.value)}
        />
        <button type="submit">Submit</button>
      </form>
      {!readData ? <div className="dots">Loading...</div> : null}
      <canvas ref={snpcovref} width={width} height={snpcovheight} />
      <canvas ref={ref} width={width} height={height} />
    </div>
  );
}

export default App;
